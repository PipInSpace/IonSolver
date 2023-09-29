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



__kernel void add(__global float* buffer, float scalar) {
    uint index = get_global_id(0) + (get_global_id(1) * get_global_size(1));
    buffer[index] += scalar;
}

__kernel void addtwo(__global float* buffer, float scalar) {
    uint index = get_global_id(0) + (get_global_id(1) * get_global_size(1));
    buffer[index] += scalar;
}

__kernel void gauss_seidel_step(__global float* buffer_new, __global float* buffer_old, float factor_a, float factor_c) {
    uint x = get_global_id(0);
    uint y = get_global_id(1);
    uint n_x = get_global_size(0);
    uint n_y = get_global_size(1);

    // Linear solving is impossible on buffer edges
    if (x != 0 && y != 0 && x != n_x-1 && y != n_y-1) {
        buffer_new[x + y*n_y] = (buffer_old[x + y*n_y] + 
        factor_a * (buffer_new[x - 1 + y*n_y] + buffer_new[x + 1 + y*n_y] + buffer_new[x + (y-1)*n_y] + buffer_new[x + (y+1)*n_y])) 
        / factor_c;
    }
}

__kernel void advect(int b, __global float* buffer_new, __global float* buffer_old, __global float* force_x, __global float* force_y, float dt) {
    uint x_index = get_global_id(0);
    uint y_index = get_global_id(1);
    uint n_x = get_global_size(0);
    uint n_y = get_global_size(1);

    // Advection is impossible at buffer edges
    // Valid x/y value indexes range from 1 to (n_x/n_y-2)
    if (x_index != 0 && y_index != 0 && x_index != n_x-1 && y_index != n_y-1) {
        uint buffer_index = x_index + (y_index * n_y);

        float x = as_float(x_index) - dt * force_x[buffer_index];
        float y = as_float(y_index) - dt * force_y[buffer_index];

        // BRANCHING VERY BAD PLS FIX
        if (x < 0.5f) {
            x = 0.5f;
        }
        if (x > as_float(n_x-2) + 0.5f) {
            x = as_float(n_x-2) + 0.5f;
        }
        uint x_idx0 = as_uint(x);
        uint x_idx1 = x_idx0 + 1;

        if (y < 0.5f) {
            y = 0.5f;
        }
        if (y > as_float(n_y-2) + 0.5f) {
            y = as_float(n_y-2) + 0.5f;
        }
        uint y_idx0 = as_uint(y);
        uint y_idx1 = y_idx0 + 1;

        float fac_s1 = x - as_float(x_idx0);
        float fac_s0 = 1.0f - fac_s1;
        float fac_t1 = y - as_float(y_idx0);
        float fac_t0 = 1.0f - fac_t1;

        buffer_new[buffer_index] = fac_s0 * (fac_t0 * buffer_old[x_idx0 + y_idx0*n_y] + fac_t1 * buffer_old[x_idx0 + y_idx1*n_y])
        + fac_s1 * (fac_t0 * buffer_old[x_idx1 + y_idx0*n_y] + fac_t1 * buffer_old[x_idx1 + y_idx1*n_y]);
    }
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