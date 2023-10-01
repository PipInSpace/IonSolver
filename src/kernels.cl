//#if defined(FP16S) || defined(FP16C)
//#define fpxx ushort
//#else // FP32
//#define fpxx float
//#endif // FP32
//
////THESE ARE TO BE REMOVED
//#define def_Nx 2u
//#define def_Ny 2u
//#define def_Nz 2u
//#define def_N 2ul
//
//#define def_Dx 1u
//#define def_Dy 1u
//#define def_Dz 1u
//
//#define def_Ox 1 // offsets are signed integer!
//#define def_Oy 1
//#define def_Oz 1
//
//#define def_Ax 1u
//#define def_Ay 1u
//#define def_Az 1u
//
//#define def_domain_offset_x 0.0f
//#define def_domain_offset_y 0.0f
//#define def_domain_offset_z 0.0f
//
//#define D "D2Q9" // D2Q9/D3Q15/D3Q19/D3Q27
//#define def_velocity_set 9u // LBM velocity set (D2Q9/D3Q15/D3Q19/D3Q27)
//#define def_dimensions 2u // number spatial dimensions (2D or 3D)
//#define def_transfers 3u // number of DDFs that are transferred between multiple domains
//
//#define def_c 0.57735027f // lattice speed of sound c = 1/sqrt(3)*dt
//#define def_w 2.0f // relaxation rate w = dt/tau = dt/(nu/c^2+dt/2) = 1/(3*nu+1/2)
//#define def_w0 (1.0f/2.25f)
//#define def_ws (1.0f/9.0f)
//#define def_we (1.0f/36.0f)
//
//#define TYPE_S 0x01 // 0b00000001 // (stationary or moving) solid boundary
//#define TYPE_E 0x02 // 0b00000010 // equilibrium boundary (inflow/outflow)
//#define TYPE_T 0x04 // 0b00000100 // temperature boundary
//#define TYPE_F 0x08 // 0b00001000 // fluid
//#define TYPE_I 0x10 // 0b00010000 // interface
//#define TYPE_G 0x20 // 0b00100000 // gas
//#define TYPE_X 0x40 // 0b01000000 // reserved type X
//#define TYPE_Y 0x80 // 0b10000000 // reserved type Y
//#define TYPE_MS 0x03 // 0b00000011 // cell next to moving solid boundary
//#define TYPE_BO 0x03 // 0b00000011 // any flag bit used for boundaries (temperature excluded)
//#define TYPE_IF 0x18 // 0b00011000 // change from interface to fluid
//#define TYPE_IG 0x30 // 0b00110000 // change from interface to gas
//#define TYPE_GI 0x38 // 0b00111000 // change from gas to interface
//#define TYPE_SU 0x38 // 0b00111000 // any flag bit used for SURFACE

uint3 coordinates(const uint n) { // disassemble 1D index to 3D coordinates (n -> x,y,z)
	const uint t = n%(def_Nx*def_Ny);
	return (uint3)(t%def_Nx, t/def_Nx, n/(def_Nx*def_Ny)); // n = x+(y+z*Ny)*Nx
}

float sq(const float x) {
	return x*x;
}

bool is_halo(const uint n) {
	const uint3 xyz = coordinates(n);
	return ((def_Dx>1u)&(xyz.x==0u||xyz.x>=def_Nx-1u))||((def_Dy>1u)&(xyz.y==0u||xyz.y>=def_Ny-1u))||((def_Dz>1u)&(xyz.z==0u||xyz.z>=def_Nz-1u));
}

bool is_halo_q(const uint3 xyz) {
	return ((def_Dx>1u)&(xyz.x==0u||xyz.x>=def_Nx-2u))||((def_Dy>1u)&(xyz.y==0u||xyz.y>=def_Ny-2u))||((def_Dz>1u)&(xyz.z==0u||xyz.z>=def_Nz-2u)); // halo data is kept up-to-date, so allow using halo data for rendering
}

ulong index_f(const uint n, const uint i) { // 64-bit indexing (maximum 2^32 lattice points (1624^3 lattice resolution, 225GB)
	return (ulong)i*def_N+(ulong)n; // SoA (229% faster on GPU)
}

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

void calculate_rho_u(const float* f, float* rhon, float* uxn, float* uyn, float* uzn) {
    float rho=f[0], ux, uy, uz;
    for(uint i=1u; i<def_velocity_set; i++) rho += f[i]; // calculate density from fi
    rho += 1.0f; // add 1.0f last to avoid digit extinction effects when summing up fi (perturbation method / DDF-shifting)
}

void calculate_f_eq(const float rho, float ux, float uy, float uz, float* feq) {
    const float c3=-3.0f*(sq(ux)+sq(uy)+sq(uz)), rhom1=rho-1.0f; // c3 = -2*sq(u)/(2*sq(c)), rhom1 is arithmetic optimization to minimize digit extinction
    ux *= 3.0f;
    uy *= 3.0f;
    uz *= 3.0f;
    feq[ 0] = def_w0*fma(rho, 0.5f*c3, rhom1); // 000 (identical for all velocity sets)
    //if defined d2q9
    const float u0=ux+uy, u1=ux-uy; // these pre-calculations make manual unrolling require less FLOPs
    const float rhos=def_ws*rho, rhoe=def_we*rho, rhom1s=def_ws*rhom1, rhom1e=def_we*rhom1;
    feq[ 1] = fma(rhos, fma(0.5f, fma(ux, ux, c3), ux), rhom1s); feq[ 2] = fma(rhos, fma(0.5f, fma(ux, ux, c3), -ux), rhom1s); // +00 -00
    feq[ 3] = fma(rhos, fma(0.5f, fma(uy, uy, c3), uy), rhom1s); feq[ 4] = fma(rhos, fma(0.5f, fma(uy, uy, c3), -uy), rhom1s); // 0+0 0-0
    feq[ 5] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), u0), rhom1e); feq[ 6] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), -u0), rhom1e); // ++0 --0
    feq[ 7] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), u1), rhom1e); feq[ 8] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), -u1), rhom1e); // +-0 -+0
}

void calculate_indices(const uint n, uint* x0, uint* xp, uint* xm, uint* y0, uint* yp, uint* ym, uint* z0, uint* zp, uint* zm) {
    const uint3 xyz = coordinates(n);
    *x0 =   xyz.x; // pre-calculate indices (periodic boundary conditions)
    *xp =  (xyz.x       +1u)%def_Nx;
    *xm =  (xyz.x+def_Nx-1u)%def_Nx;
    *y0 =   xyz.y                   *def_Nx;
    *yp = ((xyz.y       +1u)%def_Ny)*def_Nx;
    *ym = ((xyz.y+def_Ny-1u)%def_Ny)*def_Nx;
    *z0 =   xyz.z                   *def_Ny*def_Nx;
    *zp = ((xyz.z       +1u)%def_Nz)*def_Ny*def_Nx;
    *zm = ((xyz.z+def_Nz-1u)%def_Nz)*def_Ny*def_Nx;
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