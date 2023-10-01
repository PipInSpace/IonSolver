#if defined(FP16S) || defined(FP16C)
#define fpxx ushort
#else // FP32
#define fpxx float
#endif // FP32

//THESE ARE TO BE REMOVED
#define def_Nx 2u
#define def_Ny 2u
#define def_Nz 2u
#define def_N 2ul

#define def_Dx "+to_string(Dx)+"u
#define def_Dy "+to_string(Dy)+"u
#define def_Dz "+to_string(Dz)+"u

#define def_Ox "+to_string(Ox)+" // offsets are signed integer!
#define def_Oy "+to_string(Oy)+"
#define def_Oz "+to_string(Oz)+"

#define def_Ax "+to_string(Ny*Nz)+"u
#define def_Ay "+to_string(Nz*Nx)+"u
#define def_Az "+to_string(Nx*Ny)+"u

#define def_domain_offset_x "+to_string((float)Ox+(float)(Dx>1u)-0.5f*((float)Dx-1.0f)*(float)(Nx-2u*(Dx>1u)))+"f
#define def_domain_offset_y "+to_string((float)Oy+(float)(Dy>1u)-0.5f*((float)Dy-1.0f)*(float)(Ny-2u*(Dy>1u)))+"f
#define def_domain_offset_z "+to_string((float)Oz+(float)(Dz>1u)-0.5f*((float)Dz-1.0f)*(float)(Nz-2u*(Dz>1u)))+"f

#define D to_string(dimensions)+"Q"+to_string(velocity_set)+" // D2Q9/D3Q15/D3Q19/D3Q27
#define def_velocity_set 2u // LBM velocity set (D2Q9/D3Q15/D3Q19/D3Q27)
#define def_dimensions "+to_string(dimensions)+"u // number spatial dimensions (2D or 3D)
#define def_transfers "+to_string(transfers)+"u // number of DDFs that are transferred between multiple domains

#define def_c 0.57735027f // lattice speed of sound c = 1/sqrt(3)*dt
#define def_w 2.0f // relaxation rate w = dt/tau = dt/(nu/c^2+dt/2) = 1/(3*nu+1/2)

#define TYPE_S 0x01 // 0b00000001 // (stationary or moving) solid boundary
#define TYPE_E 0x02 // 0b00000010 // equilibrium boundary (inflow/outflow)
#define TYPE_T 0x04 // 0b00000100 // temperature boundary
#define TYPE_F 0x08 // 0b00001000 // fluid
#define TYPE_I 0x10 // 0b00010000 // interface
#define TYPE_G 0x20 // 0b00100000 // gas
#define TYPE_X 0x40 // 0b01000000 // reserved type X
#define TYPE_Y 0x80 // 0b10000000 // reserved type Y
#define TYPE_MS 0x03 // 0b00000011 // cell next to moving solid boundary
#define TYPE_BO 0x03 // 0b00000011 // any flag bit used for boundaries (temperature excluded)
#define TYPE_IF 0x18 // 0b00011000 // change from interface to fluid
#define TYPE_IG 0x30 // 0b00110000 // change from interface to gas
#define TYPE_GI 0x38 // 0b00111000 // change from gas to interface
#define TYPE_SU 0x38 // 0b00111000 // any flag bit used for SURFACE

void store_f(const uint n, const float* fhn, global fpxx* fi, const uint* j, const ulong t) {
	store(fi, index_f(n, 0u), fhn[0]); // Esoteric-Pull
	for(uint i=1u; i<def_velocity_set; i+=2u) {
		store(fi, index_f(j[i], t%2ul ? i+1u : i   ), fhn[i   ]);
		store(fi, index_f(n   , t%2ul ? i    : i+1u), fhn[i+1u]);
    }
}

void load_f(const uint n, float* fhn, const global fpxx* fi, const uint* j, const ulong t) {
	fhn[0] = load(fi, index_f(n, 0u)); // Esoteric-Pull
	for(uint i=1u; i<def_velocity_set; i+=2u) {
		fhn[i   ] = load(fi, index_f(n   , t%2ul ? i    : i+1u));
		fhn[i+1u] = load(fi, index_f(j[i], t%2ul ? i+1u : i   ));
	}
}

void calculate_rho_u() {
    
}

void neighbors(const uint n, uint* j) {
    uint x0, xp, xm, y0, yp, ym, z0, zp, zm;
    calculate_indices(n, &x0, &xp, &xm, &y0, &yp, &ym, &z0, &zp, &zm);
    j[0] = n;
    //TODO: if defined d2q9
    j[ 1] = xp+y0; j[ 2] = xm+y0; // +00 -00
    j[ 3] = x0+yp; j[ 4] = x0+ym; // 0+0 0-0
    j[ 5] = xp+yp; j[ 6] = xm+ym; // ++0 --0
    j[ 7] = xp+ym; j[ 8] = xm+yp; // +-0 -+0
}

__kernel void stream_collide(global fpxx* fi, global float* rho, global float* u, global uchar* flags, const ulong t, const float fx, const float fy, const float fz) {
    const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
    if(n>=(uint)def_N||is_halo(n)) return; // don't execute stream_collide() on halo
    const uchar flagsn = flags[n]; // cache flags[n] for multiple readings
    const uchar flagsn_bo=flagsn&TYPE_BO, flagsn_su=flagsn&TYPE_SU; // extract boundary and surface flags
    if(flagsn_bo==TYPE_S||flagsn_su==TYPE_G) return; // if cell is solid boundary or gas, just return

    uint j[def_velocity_set]; // neighbor indices
    neighbors(n, j); // calculate neighbor indices

    float fhn[def_velocity_set]; // local DDFs
    load_f(n, fhn, fi, j, t); // perform streaming (part 2)

    float rhon, uxn, uyn, uzn; // calculate local density and velocity for collision
    calculate_rho_u(fhn, &rhon, &uxn, &uyn, &uzn); // calculate density and velocity fields from fi

    float fxn=fx, fyn=fy, fzn=fz; // force starts as constant volume force, can be modified before call of calculate_forcing_terms(...)
    float Fin[def_velocity_set]; // forcing terms

    uxn = clamp(uxn, -def_c, def_c); // limit velocity (for stability purposes)
    uyn = clamp(uyn, -def_c, def_c); // force term: F*dt/(2*rho)
    uzn = clamp(uzn, -def_c, def_c);
    for(uint i=0u; i<def_velocity_set; i++) Fin[i] = 0.0f;

    float feq[def_velocity_set]; // equilibrium DDFs
    calculate_f_eq(rhon, uxn, uyn, uzn, feq); // calculate equilibrium DDFs
    float w = def_w; // LBM relaxation rate w = dt/tau = dt/(nu/c^2+dt/2) = 1/(3*nu+1/2)

    #if defined(SRT)
        for(uint i=0u; i<def_velocity_set; i++) fhn[i] = fma(1.0f-w, fhn[i], fma(w, feq[i], Fin[i])); // perform collision (SRT)
    #endif

    store_f(n, fhn, fi, j, t); // perform streaming (part 1)
} // stream_collide()

__kernel void initialize(global fpxx* fi, const global float* rho, global float* u, global uchar* flags) {
    const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
    if(n>=(uint)def_N||is_halo(n)) return; // don't execute initialize() on halo
    uchar flagsn = flags[n];
    const uchar flagsn_bo = flagsn&TYPE_BO; // extract boundary flags
    uint j[def_velocity_set]; // neighbor indices
    neighbors(n, j); // calculate neighbor indices
    uchar flagsj[def_velocity_set]; // cache neighbor flags for multiple readings
    for(uint i=1u; i<def_velocity_set; i++) flagsj[i] = flags[j[i]];
    if(flagsn_bo==TYPE_S) { // cell is solid
	    bool TYPE_ONLY_S = true; // has only solid neighbors
	    for(uint i=1u; i<def_velocity_set; i++) TYPE_ONLY_S = TYPE_ONLY_S&&(flagsj[i]&TYPE_BO)==TYPE_S;
	    if(TYPE_ONLY_S) {
	    	u[                 n] = 0.0f; // reset velocity for solid lattice points with only boundary neighbors
	    	u[    def_N+(ulong)n] = 0.0f;
	    	u[2ul*def_N+(ulong)n] = 0.0f;
	    }
        if(flagsn_bo==TYPE_S) {
	        u[                 n] = 0.0f; // reset velocity for all solid lattice points
	        u[    def_N+(ulong)n] = 0.0f;
	        u[2ul*def_N+(ulong)n] = 0.0f;
        }
    }
    float feq[def_velocity_set]; // f_equilibrium
    calculate_f_eq(rho[n], u[n], u[def_N+(ulong)n], u[2ul*def_N+(ulong)n], feq);
    store_f(n, feq, fi, j, 1ul); // write to fi
} // initialize()

kernel void update_fields(const global fpxx* fi, global float* rho, global float* u, const global uchar* flags, const ulong t, const float fx, const float fy, const float fz) {
    const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
    if(n>=(uint)def_N||is_halo(n)) return; // don't execute update_fields() on halo
    const uchar flagsn = flags[n];
    const uchar flagsn_bo=flagsn&TYPE_BO, flagsn_su=flagsn&TYPE_SU; // extract boundary and surface flags
    if(flagsn_bo==TYPE_S||flagsn_su==TYPE_G) return; // don't update fields for boundary or gas lattice points

    uint j[def_velocity_set]; // neighbor indices
    neighbors(n, j); // calculate neighbor indices
    float fhn[def_velocity_set]; // local DDFs
    load_f(n, fhn, fi, j, t); // perform streaming (part 2)

    float rhon, uxn, uyn, uzn; // calculate local density and velocity for collision
    calculate_rho_u(fhn, &rhon, &uxn, &uyn, &uzn); // calculate density and velocity fields from fi
    float fxn=fx, fyn=fy, fzn=fz; // force starts as constant volume force, can be modified before call of calculate_forcing_terms(...)
    {
        uxn = clamp(uxn, -def_c, def_c); // limit velocity (for stability purposes)
        uyn = clamp(uyn, -def_c, def_c); // force term: F*dt/(2*rho)
        uzn = clamp(uzn, -def_c, def_c);
    }

    rho[               n] = rhon; // update density field
    u[                 n] = uxn; // update velocity field
    u[    def_N+(ulong)n] = uyn;
    u[2ul*def_N+(ulong)n] = uzn;


}