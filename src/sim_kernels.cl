//THESE ARE TO BE REMOVED
#if defined(FP16S) || defined(FP16C)
#define fpxx ushort
#else // FP32
#define fpxx float
#endif // FP32

#define def_Nx 20u
#define def_Ny 20u
#define def_Nz 20u
#define def_N 8000ul

#define def_Dx 1u
#define def_Dy 1u
#define def_Dz 1u

#define def_Ax 1u
#define def_Ay 1u
#define def_Az 1u

#define D "D2Q9" // D2Q9/D3Q15/D3Q19/D3Q27
#define def_velocity_set 9u // LBM velocity set (D2Q9/D3Q15/D3Q19/D3Q27)
#define def_dimensions 2u // number spatial dimensions (2D or 3D)
#define def_transfers 3u // number of DDFs that are transferred between multiple domains

#define def_c 0.57735027f // lattice speed of sound c = 1/sqrt(3)*dt
#define def_w 2.0f // relaxation rate w = dt/tau = dt/(nu/c^2+dt/2) = 1/(3*nu+1/2)
#define def_w0 (1.0f/2.25f)
#define def_ws (1.0f/9.0f)
#define def_we (1.0f/36.0f)
#define def_ke 8.9875517923E9f
#define def_kmu 0.0f
#define def_charge 0.1f // Electric charge of a cell
#define def_ind_r 5 // Range of induction fill around cell
#define def_w_Q 0.1f

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

#define fpxx_copy ushort
#define load(p,o) half_to_float_custom(p[o])
#define store(p,o,x) p[o]=float_to_half_custom(x)

#define EQUILIBRIUM_BOUNDARIES
#define VOLUME_FORCE
#define MAGNETO_HYDRO
//These defines are for code completion only and are removed from the code before compilation 
#define EndTempDefines%

// Helper functions 
float sq(const float x) {
	return x*x;
}
uint uint_sq(const uint x) {
	return x*x;
}
float cb(const float x) {
	return x*x*x;
}
ushort float_to_half_custom(const float x) { // custom 16-bit floating-point format, 1-4-11, exp-15, +-1.99951168, +-6.10351562E-5, +-2.98023224E-8, 3.612 digits
	const uint b = as_uint(x)+0x00000800; // round-to-nearest-even: add last bit after truncated mantissa
	const uint e = (b&0x7F800000)>>23; // exponent
	const uint m = b&0x007FFFFF; // mantissa; in line below: 0x007FF800 = 0x00800000-0x00000800 = decimal indicator flag - initial rounding
	return (b&0x80000000)>>16 | (e>112)*((((e-112)<<11)&0x7800)|m>>12) | ((e<113)&(e>100))*((((0x007FF800+m)>>(124-e))+1)>>1); // sign : normalized : denormalized (assume [-2,2])
}
float half_to_float_custom(const ushort x) { // custom 16-bit floating-point format, 1-4-11, exp-15, +-1.99951168, +-6.10351562E-5, +-2.98023224E-8, 3.612 digits
	const uint e = (x&0x7800)>>11; // exponent
	const uint m = (x&0x07FF)<<12; // mantissa
	const uint v = as_uint((float)m)>>23; // evil log2 bit hack to count leading zeros in denormalized format
	return as_float((x&0x8000)<<16 | (e!=0)*((e+112)<<23|m) | ((e==0)&(m!=0))*((v-37)<<23|((m<<(150-v))&0x007FF000))); // sign : normalized : denormalized
}
// cube of magnitude of v
float sqmagnitude(uint3 v){
	return sq(v.x) + sq(v.y) + sq(v.z);
}
float cbmagnitude(uint3 v){
	return cb(sqrt((float)(uint_sq(v.x) + uint_sq(v.y) + uint_sq(v.z))));
}
int int_max(int x, int y) {
	if (x > y) {
		return x;
	} else {
		return y;
	}
}
int int_min(int x, int y) {
	if (x < y) {
		return x;
	} else {
		return y;
	}
}


float3 position(const uint3 xyz) { // 3D coordinates to 3D position
	return (float3)((float)xyz.x+0.5f-0.5f*(float)def_Nx, (float)xyz.y+0.5f-0.5f*(float)def_Ny, (float)xyz.z+0.5f-0.5f*(float)def_Nz);
}
uint3 coordinates(const uint n) { // disassemble 1D index to 3D coordinates (n -> x,y,z)
	const uint t = n%(def_Nx*def_Ny);
	return (uint3)(t%def_Nx, t/def_Nx, n/(def_Nx*def_Ny)); // n = x+(y+z*Ny)*Nx
}
uint index(const uint3 xyz) { // assemble 1D index from 3D coordinates (x,y,z -> n)
	return xyz.x+(xyz.y+xyz.z*def_Ny)*def_Nx; // n = x+(y+z*Ny)*Nx
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
void calculate_f_eq(const float rho, float ux, float uy, float uz, float* feq) {
    const float c3=-3.0f*(sq(ux)+sq(uy)+sq(uz)), rhom1=rho-1.0f; // c3 = -2*sq(u)/(2*sq(c)), rhom1 is arithmetic optimization to minimize digit extinction
    ux *= 3.0f;
    uy *= 3.0f;
    uz *= 3.0f;
    feq[ 0] = def_w0*fma(rho, 0.5f*c3, rhom1); // 000 (identical for all velocity sets)
    #if defined(D2Q9)
    const float u0=ux+uy, u1=ux-uy; // these pre-calculations make manual unrolling require less FLOPs
    const float rhos=def_ws*rho, rhoe=def_we*rho, rhom1s=def_ws*rhom1, rhom1e=def_we*rhom1;
    feq[ 1] = fma(rhos, fma(0.5f, fma(ux, ux, c3), ux), rhom1s); feq[ 2] = fma(rhos, fma(0.5f, fma(ux, ux, c3), -ux), rhom1s); // +00 -00
    feq[ 3] = fma(rhos, fma(0.5f, fma(uy, uy, c3), uy), rhom1s); feq[ 4] = fma(rhos, fma(0.5f, fma(uy, uy, c3), -uy), rhom1s); // 0+0 0-0
    feq[ 5] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), u0), rhom1e); feq[ 6] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), -u0), rhom1e); // ++0 --0
    feq[ 7] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), u1), rhom1e); feq[ 8] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), -u1), rhom1e); // +-0 -+0
    #elif defined(D3Q15)
    const float u0=ux+uy+uz, u1=ux+uy-uz, u2=ux-uy+uz, u3=-ux+uy+uz;
    const float rhos=def_ws*rho, rhoc=def_wc*rho, rhom1s=def_ws*rhom1, rhom1c=def_wc*rhom1;
    feq[ 1] = fma(rhos, fma(0.5f, fma(ux, ux, c3), ux), rhom1s); feq[ 2] = fma(rhos, fma(0.5f, fma(ux, ux, c3), -ux), rhom1s); // +00 -00
    feq[ 3] = fma(rhos, fma(0.5f, fma(uy, uy, c3), uy), rhom1s); feq[ 4] = fma(rhos, fma(0.5f, fma(uy, uy, c3), -uy), rhom1s); // 0+0 0-0
    feq[ 5] = fma(rhos, fma(0.5f, fma(uz, uz, c3), uz), rhom1s); feq[ 6] = fma(rhos, fma(0.5f, fma(uz, uz, c3), -uz), rhom1s); // 00+ 00-
    feq[ 7] = fma(rhoc, fma(0.5f, fma(u0, u0, c3), u0), rhom1c); feq[ 8] = fma(rhoc, fma(0.5f, fma(u0, u0, c3), -u0), rhom1c); // +++ ---
    feq[ 9] = fma(rhoc, fma(0.5f, fma(u1, u1, c3), u1), rhom1c); feq[10] = fma(rhoc, fma(0.5f, fma(u1, u1, c3), -u1), rhom1c); // ++- --+
    feq[11] = fma(rhoc, fma(0.5f, fma(u2, u2, c3), u2), rhom1c); feq[12] = fma(rhoc, fma(0.5f, fma(u2, u2, c3), -u2), rhom1c); // +-+ -+-
    feq[13] = fma(rhoc, fma(0.5f, fma(u3, u3, c3), u3), rhom1c); feq[14] = fma(rhoc, fma(0.5f, fma(u3, u3, c3), -u3), rhom1c); // -++ +--
    #elif defined(D3Q19)
    const float u0=ux+uy, u1=ux+uz, u2=uy+uz, u3=ux-uy, u4=ux-uz, u5=uy-uz;
    const float rhos=def_ws*rho, rhoe=def_we*rho, rhom1s=def_ws*rhom1, rhom1e=def_we*rhom1;
    feq[ 1] = fma(rhos, fma(0.5f, fma(ux, ux, c3), ux), rhom1s); feq[ 2] = fma(rhos, fma(0.5f, fma(ux, ux, c3), -ux), rhom1s); // +00 -00
    feq[ 3] = fma(rhos, fma(0.5f, fma(uy, uy, c3), uy), rhom1s); feq[ 4] = fma(rhos, fma(0.5f, fma(uy, uy, c3), -uy), rhom1s); // 0+0 0-0
    feq[ 5] = fma(rhos, fma(0.5f, fma(uz, uz, c3), uz), rhom1s); feq[ 6] = fma(rhos, fma(0.5f, fma(uz, uz, c3), -uz), rhom1s); // 00+ 00-
    feq[ 7] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), u0), rhom1e); feq[ 8] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), -u0), rhom1e); // ++0 --0
    feq[ 9] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), u1), rhom1e); feq[10] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), -u1), rhom1e); // +0+ -0-
    feq[11] = fma(rhoe, fma(0.5f, fma(u2, u2, c3), u2), rhom1e); feq[12] = fma(rhoe, fma(0.5f, fma(u2, u2, c3), -u2), rhom1e); // 0++ 0--
    feq[13] = fma(rhoe, fma(0.5f, fma(u3, u3, c3), u3), rhom1e); feq[14] = fma(rhoe, fma(0.5f, fma(u3, u3, c3), -u3), rhom1e); // +-0 -+0
    feq[15] = fma(rhoe, fma(0.5f, fma(u4, u4, c3), u4), rhom1e); feq[16] = fma(rhoe, fma(0.5f, fma(u4, u4, c3), -u4), rhom1e); // +0- -0+
    feq[17] = fma(rhoe, fma(0.5f, fma(u5, u5, c3), u5), rhom1e); feq[18] = fma(rhoe, fma(0.5f, fma(u5, u5, c3), -u5), rhom1e); // 0+- 0-+
    #elif defined(D3Q27)
    const float u0=ux+uy, u1=ux+uz, u2=uy+uz, u3=ux-uy, u4=ux-uz, u5=uy-uz, u6=ux+uy+uz, u7=ux+uy-uz, u8=ux-uy+uz, u9=-ux+uy+uz;
    const float rhos=def_ws*rho, rhoe=def_we*rho, rhoc=def_wc*rho, rhom1s=def_ws*rhom1, rhom1e=def_we*rhom1, rhom1c=def_wc*rhom1;
    feq[ 1] = fma(rhos, fma(0.5f, fma(ux, ux, c3), ux), rhom1s); feq[ 2] = fma(rhos, fma(0.5f, fma(ux, ux, c3), -ux), rhom1s); // +00 -00
    feq[ 3] = fma(rhos, fma(0.5f, fma(uy, uy, c3), uy), rhom1s); feq[ 4] = fma(rhos, fma(0.5f, fma(uy, uy, c3), -uy), rhom1s); // 0+0 0-0
    feq[ 5] = fma(rhos, fma(0.5f, fma(uz, uz, c3), uz), rhom1s); feq[ 6] = fma(rhos, fma(0.5f, fma(uz, uz, c3), -uz), rhom1s); // 00+ 00-
    feq[ 7] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), u0), rhom1e); feq[ 8] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), -u0), rhom1e); // ++0 --0
    feq[ 9] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), u1), rhom1e); feq[10] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), -u1), rhom1e); // +0+ -0-
    feq[11] = fma(rhoe, fma(0.5f, fma(u2, u2, c3), u2), rhom1e); feq[12] = fma(rhoe, fma(0.5f, fma(u2, u2, c3), -u2), rhom1e); // 0++ 0--
    feq[13] = fma(rhoe, fma(0.5f, fma(u3, u3, c3), u3), rhom1e); feq[14] = fma(rhoe, fma(0.5f, fma(u3, u3, c3), -u3), rhom1e); // +-0 -+0
    feq[15] = fma(rhoe, fma(0.5f, fma(u4, u4, c3), u4), rhom1e); feq[16] = fma(rhoe, fma(0.5f, fma(u4, u4, c3), -u4), rhom1e); // +0- -0+
    feq[17] = fma(rhoe, fma(0.5f, fma(u5, u5, c3), u5), rhom1e); feq[18] = fma(rhoe, fma(0.5f, fma(u5, u5, c3), -u5), rhom1e); // 0+- 0-+
    feq[19] = fma(rhoc, fma(0.5f, fma(u6, u6, c3), u6), rhom1c); feq[20] = fma(rhoc, fma(0.5f, fma(u6, u6, c3), -u6), rhom1c); // +++ ---
    feq[21] = fma(rhoc, fma(0.5f, fma(u7, u7, c3), u7), rhom1c); feq[22] = fma(rhoc, fma(0.5f, fma(u7, u7, c3), -u7), rhom1c); // ++- --+
    feq[23] = fma(rhoc, fma(0.5f, fma(u8, u8, c3), u8), rhom1c); feq[24] = fma(rhoc, fma(0.5f, fma(u8, u8, c3), -u8), rhom1c); // +-+ -+-
    feq[25] = fma(rhoc, fma(0.5f, fma(u9, u9, c3), u9), rhom1c); feq[26] = fma(rhoc, fma(0.5f, fma(u9, u9, c3), -u9), rhom1c); // -++ +--
    #endif
}
void calculate_rho_u(const float* f, float* rhon, float* uxn, float* uyn, float* uzn) {
    float rho=f[0], ux, uy, uz;
    for(uint i=1u; i<def_velocity_set; i++) rho += f[i]; // calculate density from fi
    rho += 1.0f; // add 1.0f last to avoid digit extinction effects when summing up fi (perturbation method / DDF-shifting)
    #if defined(D2Q9)
    ux = f[1]-f[2]+f[5]-f[6]+f[7]-f[8]; // calculate velocity from fi (alternating + and - for best accuracy)
    uy = f[3]-f[4]+f[5]-f[6]+f[8]-f[7];
    uz = 0.0f;
    #elif defined(D3Q15)
    ux = f[ 1]-f[ 2]+f[ 7]-f[ 8]+f[ 9]-f[10]+f[11]-f[12]+f[14]-f[13]; // calculate velocity from fi (alternating + and - for best accuracy)
    uy = f[ 3]-f[ 4]+f[ 7]-f[ 8]+f[ 9]-f[10]+f[12]-f[11]+f[13]-f[14];
    uz = f[ 5]-f[ 6]+f[ 7]-f[ 8]+f[10]-f[ 9]+f[11]-f[12]+f[13]-f[14];
    #elif defined(D3Q19)
    ux = f[ 1]-f[ 2]+f[ 7]-f[ 8]+f[ 9]-f[10]+f[13]-f[14]+f[15]-f[16]; // calculate velocity from fi (alternating + and - for best accuracy)
    uy = f[ 3]-f[ 4]+f[ 7]-f[ 8]+f[11]-f[12]+f[14]-f[13]+f[17]-f[18];
    uz = f[ 5]-f[ 6]+f[ 9]-f[10]+f[11]-f[12]+f[16]-f[15]+f[18]-f[17];
    #elif defined(D3Q27)
    ux = f[ 1]-f[ 2]+f[ 7]-f[ 8]+f[ 9]-f[10]+f[13]-f[14]+f[15]-f[16]+f[19]-f[20]+f[21]-f[22]+f[23]-f[24]+f[26]-f[25]; // calculate velocity from fi (alternating + and - for best accuracy)
    uy = f[ 3]-f[ 4]+f[ 7]-f[ 8]+f[11]-f[12]+f[14]-f[13]+f[17]-f[18]+f[19]-f[20]+f[21]-f[22]+f[24]-f[23]+f[25]-f[26];
    uz = f[ 5]-f[ 6]+f[ 9]-f[10]+f[11]-f[12]+f[16]-f[15]+f[18]-f[17]+f[19]-f[20]+f[22]-f[21]+f[23]-f[24]+f[25]-f[26];
    #endif
    *rhon = rho;
    *uxn = ux/rho;
    *uyn = uy/rho;
    *uzn = uz/rho;
} // calculate_rho_u
void load_f(const uint n, float* fhn, const global fpxx* fi, const uint* j, const ulong t) {
	fhn[0] = load(fi, index_f(n, 0u)); // Esoteric-Pull
	for(uint i=1u; i<def_velocity_set; i+=2u) {
		fhn[i   ] = load(fi, index_f(n   , t%2ul ? i    : i+1u));
		fhn[i+1u] = load(fi, index_f(j[i], t%2ul ? i+1u : i   ));
	}
}
void store_f(const uint n, const float* fhn, global fpxx* fi, const uint* j, const ulong t) {
	store(fi, index_f(n, 0u), fhn[0]); // Esoteric-Pull
	for(uint i=1u; i<def_velocity_set; i+=2u) {
		store(fi, index_f(j[i], t%2ul ? i+1u : i   ), fhn[i   ]);
		store(fi, index_f(n   , t%2ul ? i    : i+1u), fhn[i+1u]);
    }
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
    #if defined(D2Q9)
    j[ 1] = xp+y0; j[ 2] = xm+y0; // +00 -00
    j[ 3] = x0+yp; j[ 4] = x0+ym; // 0+0 0-0
    j[ 5] = xp+yp; j[ 6] = xm+ym; // ++0 --0
    j[ 7] = xp+ym; j[ 8] = xm+yp; // +-0 -+0
    #elif defined(D3Q15)
    j[ 1] = xp+y0+z0; j[ 2] = xm+y0+z0; // +00 -00
    j[ 3] = x0+yp+z0; j[ 4] = x0+ym+z0; // 0+0 0-0
    j[ 5] = x0+y0+zp; j[ 6] = x0+y0+zm; // 00+ 00-
    j[ 7] = xp+yp+zp; j[ 8] = xm+ym+zm; // +++ ---
    j[ 9] = xp+yp+zm; j[10] = xm+ym+zp; // ++- --+
    j[11] = xp+ym+zp; j[12] = xm+yp+zm; // +-+ -+-
    j[13] = xm+yp+zp; j[14] = xp+ym+zm; // -++ +--
    #elif defined(D3Q19)
    j[ 1] = xp+y0+z0; j[ 2] = xm+y0+z0; // +00 -00
    j[ 3] = x0+yp+z0; j[ 4] = x0+ym+z0; // 0+0 0-0
    j[ 5] = x0+y0+zp; j[ 6] = x0+y0+zm; // 00+ 00-
    j[ 7] = xp+yp+z0; j[ 8] = xm+ym+z0; // ++0 --0
    j[ 9] = xp+y0+zp; j[10] = xm+y0+zm; // +0+ -0-
    j[11] = x0+yp+zp; j[12] = x0+ym+zm; // 0++ 0--
    j[13] = xp+ym+z0; j[14] = xm+yp+z0; // +-0 -+0
    j[15] = xp+y0+zm; j[16] = xm+y0+zp; // +0- -0+
    j[17] = x0+yp+zm; j[18] = x0+ym+zp; // 0+- 0-+
    #elif defined(D3Q27)
    j[ 1] = xp+y0+z0; j[ 2] = xm+y0+z0; // +00 -00
    j[ 3] = x0+yp+z0; j[ 4] = x0+ym+z0; // 0+0 0-0
    j[ 5] = x0+y0+zp; j[ 6] = x0+y0+zm; // 00+ 00-
    j[ 7] = xp+yp+z0; j[ 8] = xm+ym+z0; // ++0 --0
    j[ 9] = xp+y0+zp; j[10] = xm+y0+zm; // +0+ -0-
    j[11] = x0+yp+zp; j[12] = x0+ym+zm; // 0++ 0--
    j[13] = xp+ym+z0; j[14] = xm+yp+z0; // +-0 -+0
    j[15] = xp+y0+zm; j[16] = xm+y0+zp; // +0- -0+
    j[17] = x0+yp+zm; j[18] = x0+ym+zp; // 0+- 0-+
    j[19] = xp+yp+zp; j[20] = xm+ym+zm; // +++ ---
    j[21] = xp+yp+zm; j[22] = xm+ym+zp; // ++- --+
    j[23] = xp+ym+zp; j[24] = xm+yp+zm; // +-+ -+-
    j[25] = xm+yp+zp; j[26] = xp+ym+zm; // -++ +--
    #endif
} //neighbors
float3 load_u(const uint n, const global float* u) {
	return (float3)(u[n], u[def_N+(ulong)n], u[2ul*def_N+(ulong)n]);
}
float calculate_Q_cached(const float3 u0, const float3 u1, const float3 u2, const float3 u3, const float3 u4, const float3 u5) { // Q-criterion
	const float duxdx=u0.x-u1.x, duydx=u0.y-u1.y, duzdx=u0.z-u1.z; // du/dx = (u2-u0)/2
	const float duxdy=u2.x-u3.x, duydy=u2.y-u3.y, duzdy=u2.z-u3.z;
	const float duxdz=u4.x-u5.x, duydz=u4.y-u5.y, duzdz=u4.z-u5.z;
	const float omega_xy=duxdy-duydx, omega_xz=duxdz-duzdx, omega_yz=duydz-duzdy; // antisymmetric tensor, omega_xx = omega_yy = omega_zz = 0
	const float s_xx2=duxdx, s_yy2=duydy, s_zz2=duzdz; // s_xx2 = s_xx/2, s_yy2 = s_yy/2, s_zz2 = s_zz/2
	const float s_xy=duxdy+duydx, s_xz=duxdz+duzdx, s_yz=duydz+duzdy; // symmetric tensor
	const float omega2 = sq(omega_xy)+sq(omega_xz)+sq(omega_yz); // ||omega||_2^2
	const float s2 = 2.0f*(sq(s_xx2)+sq(s_yy2)+sq(s_zz2))+sq(s_xy)+sq(s_xz)+sq(s_yz); // ||s||_2^2
	return 0.25f*(omega2-s2); // Q = 1/2*(||omega||_2^2-||s||_2^2), addidional factor 1/2 from cental finite differences of velocity
} // calculate_Q_cached()
float calculate_Q(const uint n, const global float* u) { // Q-criterion
	uint x0, xp, xm, y0, yp, ym, z0, zp, zm;
	calculate_indices(n, &x0, &xp, &xm, &y0, &yp, &ym, &z0, &zp, &zm);
	uint j[6];
	j[0] = xp+y0+z0; j[1] = xm+y0+z0; // +00 -00
	j[2] = x0+yp+z0; j[3] = x0+ym+z0; // 0+0 0-0
	j[4] = x0+y0+zp; j[5] = x0+y0+zm; // 00+ 00-
	return calculate_Q_cached(load_u(j[0], u), load_u(j[1], u), load_u(j[2], u), load_u(j[3], u), load_u(j[4], u), load_u(j[5], u));
} // calculate_Q()
float c(const uint i) { // avoid constant keyword by encapsulating data in function which gets inlined by compiler
	const float c[3u*def_velocity_set] = {
	#if defined(D2Q9)
		0, 1,-1, 0, 0, 1,-1, 1,-1, // x
		0, 0, 0, 1,-1, 1,-1,-1, 1, // y
		0, 0, 0, 0, 0, 0, 0, 0, 0  // z
	#elif defined(D3Q15)
		0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, // x
		0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1, // y
		0, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1, 1,-1  // z
	#elif defined(D3Q19)
		0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, // x
		0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, // y
		0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1  // z
	#elif defined(D3Q27)
		0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, // x
		0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, 1,-1, // y
		0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1  // z
	#endif // D3Q27
	};
	return c[i];
}
float w(const uint i) { // avoid constant keyword by encapsulating data in function which gets inlined by compiler
	const float w[def_velocity_set] = { def_w0, // velocity set weights
	#if defined(D2Q9)
		def_ws, def_ws, def_ws, def_ws, def_we, def_we, def_we, def_we
	#elif defined(D3Q15)
		def_ws, def_ws, def_ws, def_ws, def_ws, def_ws,
		def_wc, def_wc, def_wc, def_wc, def_wc, def_wc, def_wc, def_wc
	#elif defined(D3Q19)
		def_ws, def_ws, def_ws, def_ws, def_ws, def_ws,
		def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we
	#elif defined(D3Q27)
		def_ws, def_ws, def_ws, def_ws, def_ws, def_ws,
		def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we,
		def_wc, def_wc, def_wc, def_wc, def_wc, def_wc, def_wc, def_wc
	#endif
	};
	return w[i];
}
#ifdef VOLUME_FORCE
void calculate_forcing_terms(const float ux, const float uy, const float uz, const float fx, const float fy, const float fz, float* Fin) { // calculate volume force terms Fin from velocity field (Guo forcing, Krueger p.233f)
	#ifdef D2Q9
		const float uF = -0.33333334f*fma(ux, fx, uy*fy); // 2D
	#else
		const float uF = -0.33333334f*fma(ux, fx, fma(uy, fy, uz*fz)); // 3D
	#endif
	Fin[0] = 9.0f*def_w0*uF ; // 000 (identical for all velocity sets)
	for(uint i=1u; i<def_velocity_set; i++) { // loop is entirely unrolled by compiler, no unnecessary FLOPs are happening
		Fin[i] = 9.0f*w(i)*fma(c(i)*fx+c(def_velocity_set+i)*fy+c(2u*def_velocity_set+i)*fz, c(i)*ux+c(def_velocity_set+i)*uy+c(2u*def_velocity_set+i)*uz+0.33333334f, uF);
	}
}
#endif // VOLUME_FORCE

#ifdef MAGNETO_HYDRO
// Charge advection
void neighbors_charge(const uint n, uint* j7) { // calculate neighbor indices
	uint x0, xp, xm, y0, yp, ym, z0, zp, zm;
	calculate_indices(n, &x0, &xp, &xm, &y0, &yp, &ym, &z0, &zp, &zm);
	j7[0] = n;
	j7[1] = xp+y0+z0; j7[2] = xm+y0+z0; // +00 -00
	j7[3] = x0+yp+z0; j7[4] = x0+ym+z0; // 0+0 0-0
	j7[5] = x0+y0+zp; j7[6] = x0+y0+zm; // 00+ 00-
}
void calculate_q_eq(const float Q, const float ux, const float uy, const float uz, float* qeq) { // calculate q_equilibrium from density and velocity field (perturbation method / DDF-shifting)
	const float wsT4=0.5f*Q, wsTm1=0.125f*(Q-1.0f); // 0.125f*Q*4.0f (straight directions in D3Q7), wsTm1 is arithmetic optimization to minimize digit extinction, lattice speed of sound is 1/2 for D3Q7 and not 1/sqrt(3)
	qeq[0] = fma(0.25f, Q, -0.25f); // 000
	qeq[1] = fma(wsT4, ux, wsTm1); qeq[2] = fma(wsT4, -ux, wsTm1); // +00 -00, source: http://dx.doi.org/10.1016/j.ijheatmasstransfer.2009.11.014
	qeq[3] = fma(wsT4, uy, wsTm1); qeq[4] = fma(wsT4, -uy, wsTm1); // 0+0 0-0
	qeq[5] = fma(wsT4, uz, wsTm1); qeq[6] = fma(wsT4, -uz, wsTm1); // 00+ 00-
}
void load_q(const uint n, float* qhn, const global fpxx* qi, const uint* j7, const ulong t) {
	qhn[0] = load(qi, index_f(n, 0u)); // Esoteric-Pull
	for(uint i=1u; i<7u; i+=2u) {
		qhn[i   ] = load(qi, index_f(n    , t%2ul ? i    : i+1u));
		qhn[i+1u] = load(qi, index_f(j7[i], t%2ul ? i+1u : i   ));
	}
}
void store_q(const uint n, const float* qhn, global fpxx* qi, const uint* j7, const ulong t) {
	store(qi, index_f(n, 0u), qhn[0]); // Esoteric-Pull
	for(uint i=1u; i<7u; i+=2u) {
		store(qi, index_f(j7[i], t%2ul ? i+1u : i   ), qhn[i   ]);
		store(qi, index_f(n    , t%2ul ? i    : i+1u), qhn[i+1u]);
	}
}
#endif

// Simulation kernel functions
__kernel void stream_collide(global fpxx* fi, global float* rho, global float* u, global uchar* flags, const ulong t, const float fx, const float fy, const float fz 
#ifdef FORCE_FIELD
, const global float* F 
#endif // FORCE_FIELD
#ifdef MAGNETO_HYDRO
, const global float* E	// static electric field
, const global float* B	// static magnetic flux
, global float* E_dyn	// dynamic electric field
, global float* B_dyn	// dynamic magnetic flux
, global float* qi		// charge ddfs
, global float* Q		// charge
#endif // MAGNETO_HYDRO
) {
    const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
    if(n>=(uint)def_N||is_halo(n)) return; // don't execute stream_collide() on halo
    const uchar flagsn = flags[n]; // cache flags[n] for multiple readings
    const uchar flagsn_bo=flagsn&TYPE_BO, flagsn_su=flagsn&TYPE_SU; // extract boundary and surface flags
    if(flagsn_bo==TYPE_S||flagsn_su==TYPE_G) return; // if cell is solid boundary or gas, just return

    uint j[def_velocity_set]; // neighbor indices
    neighbors(n, j); // calculate neighbor indices

    float fhn[def_velocity_set]; // local DDFs
    load_f(n, fhn, fi, j, t); // perform streaming (part 2)

	ulong nxi=(ulong)n, nyi=def_N+(ulong)n, nzi=2ul*def_N+(ulong)n; // n indecies for x, y and z components

    float rhon, uxn, uyn, uzn; // calculate local density and velocity for collision

    #ifndef EQUILIBRIUM_BOUNDARIES // EQUILIBRIUM_BOUNDARIES
        calculate_rho_u(fhn, &rhon, &uxn, &uyn, &uzn); // calculate density and velocity fields from fi
    #else
        if(flagsn_bo==TYPE_E) {
        	rhon = rho[n]; // apply preset velocity/density
        	uxn  = u[nxi];
        	uyn  = u[nyi];
        	uzn  = u[nzi];
        } else {
        	calculate_rho_u(fhn, &rhon, &uxn, &uyn, &uzn); // calculate density and velocity fields from fi
        }
    #endif // EQUILIBRIUM_BOUNDARIES

    float fxn=fx, fyn=fy, fzn=fz; // force starts as constant volume force, can be modified before call of calculate_forcing_terms(...)

    float Fin[def_velocity_set]; // forcing terms

	#ifdef FORCE_FIELD
	{ // separate block to avoid variable name conflicts
		fxn += F[nxi]; // apply force field
		fyn += F[nyi];
		fzn += F[nzi];
	}
	#endif

	#ifdef MAGNETO_HYDRO
	{
		// Advection of charge. Cell charge is stored in charge ddfs 'qi' for advection with 'fi'.
		uint j7[7]; // neighbors of D3Q7 subset
		neighbors_charge(n, j7);
		float qhn[7]; // read from qA and stream to gh (D3Q7 subset, periodic boundary conditions)
		load_q(n, qhn, qi, j7, t); // perform streaming (part 2)
		float Qn = 0.0f;
		for(uint i=0u; i<7u; i++) Qn += qhn[i]; // calculate charge from q
		Qn += 1.0f; // add 1.0f last to avoid digit extinction effects when summing up qi (perturbation method / DDF-shifting)
		float qeq[7]; // cache f_equilibrium[n]
		calculate_q_eq(Qn, uxn, uyn, uzn, qeq); // calculate equilibrium DDFs

		#ifdef UPDATE_FIELDS
			Q[n] = Qn; // update charge field
		#endif // UPDATE_FIELDS

		for(uint i=0u; i<7u; i++) qhn[i] = fma(1.0f-def_w_Q, qhn[i], def_w_Q*qeq[i]); // perform collision
		store_q(n, qhn, qi, j7, t); // perform streaming (part 1)

		// F = charge * (E + (U cross B))
		fxn += Qn * (E_dyn[nxi] + uyn*B_dyn[nzi] - uzn*B_dyn[nyi]); // force = charge * (electric field + magnetic field x U)
		fyn += Qn * (E_dyn[nyi] + uzn*B_dyn[nxi] - uxn*B_dyn[nzi]);
		fzn += Qn * (E_dyn[nzi] + uxn*B_dyn[nyi] - uyn*B_dyn[nxi]);

		// Clear dynamic buffers with static buffers for recomputation
		B_dyn[nxi] = B[nxi];
		B_dyn[nyi] = B[nyi];
		B_dyn[nzi] = B[nzi];
		E_dyn[nxi] = E[nxi];
		E_dyn[nyi] = E[nyi];
		E_dyn[nzi] = E[nzi];
	}
	#endif// MAGNETO_HYDRO

	#ifdef VOLUME_FORCE
		const float rho2 = 0.5f/rhon; // apply external volume force (Guo forcing, Krueger p.233f)
		uxn = clamp(fma(fxn, rho2, uxn), -def_c, def_c); // limit velocity (for stability purposes)
		uyn = clamp(fma(fyn, rho2, uyn), -def_c, def_c); // force term: F*dt/(2*rho)
		uzn = clamp(fma(fzn, rho2, uzn), -def_c, def_c);
		calculate_forcing_terms(uxn, uyn, uzn, fxn, fyn, fzn, Fin); // calculate volume force terms Fin from velocity field (Guo forcing, Krueger p.233f)
	#else // VOLUME_FORCE
    	uxn = clamp(uxn, -def_c, def_c); // limit velocity (for stability purposes)
    	uyn = clamp(uyn, -def_c, def_c); // force term: F*dt/(2*rho)
    	uzn = clamp(uzn, -def_c, def_c);
    	for(uint i=0u; i<def_velocity_set; i++) Fin[i] = 0.0f;
	#endif // VOLUME_FORCE


	#ifndef EQUILIBRIUM_BOUNDARIES
	#ifdef UPDATE_FIELDS
		rho[n] = rhon; // update density field
		u[nxi] = uxn; // update velocity field
		u[nyi] = uyn;
		u[nzi] = uzn;
	#endif // UPDATE_FIELDS
	#else // EQUILIBRIUM_BOUNDARIES
	#ifdef UPDATE_FIELDS
		if(flagsn_bo!=TYPE_E) { // only update fields for non-TYPE_E cells
			rho[n] = rhon; // update density field
			u[nxi] = uxn; // update velocity field
			u[nyi] = uyn;
			u[nzi] = uzn;
		}
	#endif // UPDATE_FIELDS
	#endif // EQUILIBRIUM_BOUNDARIES

    float feq[def_velocity_set]; // equilibrium DDFs
    calculate_f_eq(rhon, uxn, uyn, uzn, feq); // calculate equilibrium DDFs
    float w = def_w; // LBM relaxation rate w = dt/tau = dt/(nu/c^2+dt/2) = 1/(3*nu+1/2)

    #if defined(SRT) // SRT
		#ifdef VOLUME_FORCE
			const float c_tau = fma(w, -0.5f, 1.0f);
			for(uint i=0u; i<def_velocity_set; i++) Fin[i] *= c_tau;
		#endif // VOLUME_FORCE

        #ifndef EQUILIBRIUM_BOUNDARIES
            for(uint i=0u; i<def_velocity_set; i++) fhn[i] = fma(1.0f-w, fhn[i], fma(w, feq[i], Fin[i])); // perform collision (SRT)
        #else
            for(uint i=0u; i<def_velocity_set; i++) fhn[i] = flagsn_bo==TYPE_E ? feq[i] : fma(1.0f-w, fhn[i], fma(w, feq[i], Fin[i])); // perform collision (SRT)
        #endif // EQUILIBRIUM_BOUNDARIES

    #elif defined(TRT) // TRT
        const float wp = w; // TRT: inverse of "+" relaxation time
        const float wm = 1.0f/(0.1875f/(1.0f/w-0.5f)+0.5f); // TRT: inverse of "-" relaxation time wm = 1.0f/(0.1875f/(3.0f*nu)+0.5f), nu = (1.0f/w-0.5f)/3.0f;

		#ifdef VOLUME_FORCE
			const float c_taup=fma(wp, -0.25f, 0.5f), c_taum=fma(wm, -0.25f, 0.5f); // source: https://arxiv.org/pdf/1901.08766.pdf
			float Fib[def_velocity_set]; // F_bar
			Fib[0] = Fin[0];
			for(uint i=1u; i<def_velocity_set; i+=2u) {
				Fib[i   ] = Fin[i+1u];
				Fib[i+1u] = Fin[i   ];
			}
			for(uint i=0u; i<def_velocity_set; i++) Fin[i] = fma(c_taup, Fin[i]+Fib[i], c_taum*(Fin[i]-Fib[i]));
		#endif // VOLUME_FORCE

        float fhb[def_velocity_set]; // fhn in inverse directions
        float feb[def_velocity_set]; // feq in inverse directions
        fhb[0] = fhn[0];
        feb[0] = feq[0];
        for(uint i=1u; i<def_velocity_set; i+=2u) {
        	fhb[i   ] = fhn[i+1u];
        	fhb[i+1u] = fhn[i   ];
        	feb[i   ] = feq[i+1u];
        	feb[i+1u] = feq[i   ];
        }
        #ifndef EQUILIBRIUM_BOUNDARIES
            for(uint i=0u; i<def_velocity_set; i++) fhn[i] = fma(0.5f*wp, feq[i]-fhn[i]+feb[i]-fhb[i], fma(0.5f*wm, feq[i]-feb[i]-fhn[i]+fhb[i], fhn[i]+Fin[i])); // perform collision (TRT)
        #else // EQUILIBRIUM_BOUNDARIES
            for(uint i=0u; i<def_velocity_set; i++) fhn[i] = flagsn_bo==TYPE_E ? feq[i] : fma(0.5f*wp, feq[i]-fhn[i]+feb[i]-fhb[i], fma(0.5f*wm, feq[i]-feb[i]-fhn[i]+fhb[i], fhn[i]+Fin[i])); // perform collision (TRT)
        #endif // EQUILIBRIUM_BOUNDARIES
    #endif // TRT

    store_f(n, fhn, fi, j, t); // perform streaming (part 1)
} // stream_collide()

__kernel void initialize(global fpxx* fi, global float* rho, global float* u, global uchar* flags
#ifdef MAGNETO_HYDRO
, const global float* E
, const global float* B
, global float* E_dyn
, global float* B_dyn
, global float* qi
, const global float* Q
#endif // MAGNETO_HYDRO
) {
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
	#ifdef MAGNETO_HYDRO
		// Initialize charge ddfs
		float qeq[7]; // q_equilibrium
		calculate_q_eq(Qn, u[n], u[def_N+(ulong)n], u[2ul*def_N+(ulong)n], qeq);
		uint j7[7]; // neighbors of D3Q7 subset
		neighbors_charge(n, j7);
		store_q(n, qeq, qi, j7, 1ul); // write to qi. perform streaming (part 1)

		// Clear dyn with static for recomputation
		B_dyn[                 n] = B[                 n];
		B_dyn[    def_N+(ulong)n] = B[    def_N+(ulong)n];
		B_dyn[2ul*def_N+(ulong)n] = B[2ul*def_N+(ulong)n];
		E_dyn[                 n] = E[                 n];
		E_dyn[    def_N+(ulong)n] = E[    def_N+(ulong)n];
		E_dyn[2ul*def_N+(ulong)n] = E[2ul*def_N+(ulong)n];
	#endif // MAGNETO_HYDRO
} // initialize()

__kernel void update_fields(const global fpxx* fi, global float* rho, global float* u, const global uchar* flags, const ulong t, const float fx, const float fy, const float fz) {
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
	// update charge field (later)


} // update_fields()

#ifdef MAGNETO_HYDRO
__kernel void update_e_b_dynamic(global float* E_dyn, global float* B_dyn, const global float* q, const global float* u, const global uchar* flags) {
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
    if(n>=(uint)def_N||is_halo(n)) return; // don't execute update_e_b_dynamic() on halo
    const uchar flagsn = flags[n]; // cache flags[n] for multiple readings
    const uchar flagsn_bo=flagsn&TYPE_BO, flagsn_su=flagsn&TYPE_SU; // extract boundary and surface flags
    if(flagsn_bo==TYPE_S||flagsn_su==TYPE_G) return; // if cell is solid boundary or gas, just return

	const uint3 coord_n = coordinates(n); // Cell coordinate

	const uint x_upper = int_min(coord_n.x + def_ind_r + 1, def_Nx);
	const uint y_upper = int_min(coord_n.y + def_ind_r + 1, def_Ny);
	const uint z_upper = int_min(coord_n.z + def_ind_r + 1, def_Nz);

	float3 e, b;
	
	for (uint x = int_max(coord_n.x - def_ind_r, 0); x < x_upper; x++) {
		for (uint y = int_max(coord_n.y - def_ind_r, 0); y < y_upper; y++) {
			for (uint z = int_max(coord_n.z - def_ind_r, 0); z < z_upper; z++) {

				// _c vars describe surronding cells 
				const uint n_c = x + (y + z * def_Nx) * def_Ny;
				if (n == n_c) continue;
					
				const float q_c = q[n_c]; // charge of nearby cell
				const float3 v_c = {u[n_c], u[(ulong)n_c+def_N], u[(ulong)n_c+def_N*2ul]}; // velocity of nearby cell

				// precalculation for both fields
				const uint3 vec_r = coordinates(n_c) - coord_n;
				const float3 pre_field =  convert_float3(vec_r) / cbmagnitude(vec_r);

				e += def_ke * q_c * pre_field; // E imparted by nearby cell (Coulomb)
				b += def_kmu * q_c * cross(v_c, pre_field); // B imparted by nearby cell (Biot-Savart)
			}
		}
	}

	// update buffers
	E_dyn[n						] += e.x;
	E_dyn[(ulong)n+def_N    	] += e.y;
	E_dyn[(ulong)n+def_N*2ul	] += e.z;

	B_dyn[n						] += b.x;
	B_dyn[(ulong)n+def_N		] += b.y;
	B_dyn[(ulong)n+def_N*2ul	] += b.z;
} // update_e_b_dynamic()
#endif

// Inter-Domain Transfer kernels
uint get_area(const uint direction) {
	const uint A[3] = { def_Ax, def_Ay, def_Az };
	return A[direction];
}
// Return 1D index of cell to be transferred from id a and xyz direction (0, 1, 2) for pos and neg directions
uint index_extract_p(const uint a, const uint direction) {
	const uint3 coordinates[3] = { (uint3)(def_Nx-2u, a%def_Ny, a/def_Ny), (uint3)(a/def_Nz, def_Ny-2u, a%def_Nz), (uint3)(a%def_Nx, a/def_Nx, def_Nz-2u) };
	return index(coordinates[direction]);
}
uint index_extract_m(const uint a, const uint direction) {
	const uint3 coordinates[3] = { (uint3)(       1u, a%def_Ny, a/def_Ny), (uint3)(a/def_Nz,        1u, a%def_Nz), (uint3)(a%def_Nx, a/def_Nx,        1u) };
	return index(coordinates[direction]);
}
uint index_insert_p(const uint a, const uint direction) {
	const uint3 coordinates[3] = { (uint3)(def_Nx-1u, a%def_Ny, a/def_Ny), (uint3)(a/def_Nz, def_Ny-1u, a%def_Nz), (uint3)(a%def_Nx, a/def_Nx, def_Nz-1u) };
	return index(coordinates[direction]);
}
uint index_insert_m(const uint a, const uint direction) {
	const uint3 coordinates[3] = { (uint3)(       0u, a%def_Ny, a/def_Ny), (uint3)(a/def_Nz,        0u, a%def_Nz), (uint3)(a%def_Nx, a/def_Nx,        0u) };
	return index(coordinates[direction]);
}
// Returns an index for the transferred ddfs
uint index_transfer(const uint side_i) {
	const uchar index_transfer_data[2u*def_dimensions*def_transfers] = {
	#if defined(D2Q9)
		1,  5,  7, // xp
		2,  6,  8, // xm
		3,  5,  8, // yp
		4,  6,  7  // ym
	#elif defined(D3Q15)
		1,  7, 14,  9, 11, // xp
		2,  8, 13, 10, 12, // xm
		3,  7, 12,  9, 13, // yp
		4,  8, 11, 10, 14, // ym
		5,  7, 10, 11, 13, // zp
		6,  8,  9, 12, 14  // zm
	#elif defined(D3Q19)
		1,  7, 13,  9, 15, // xp
		2,  8, 14, 10, 16, // xm
		3,  7, 14, 11, 17, // yp
		4,  8, 13, 12, 18, // ym
		5,  9, 16, 11, 18, // zp
		6, 10, 15, 12, 17  // zm
	#elif defined(D3Q27)
		1,  7, 13,  9, 15, 19, 26, 21, 23, // xp
		2,  8, 14, 10, 16, 20, 25, 22, 24, // xm
		3,  7, 14, 11, 17, 19, 24, 21, 25, // yp
		4,  8, 13, 12, 18, 20, 23, 22, 26, // ym
		5,  9, 16, 11, 18, 19, 22, 23, 25, // zp
		6, 10, 15, 12, 17, 20, 21, 24, 26  // zm
	#endif // D3Q27
	};
	return (uint)index_transfer_data[side_i];
}

// Fi
void extract_fi(const uint a, const uint A, const uint n, const uint side, const ulong t, global fpxx_copy* transfer_buffer, const global fpxx_copy* fi) {
	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	for(uint b=0u; b<def_transfers; b++) {
		const uint i = index_transfer(side*def_transfers+b);
		const ulong index = index_f(i%2u ? j[i] : n, t%2ul ? (i%2u ? i+1u : i-1u) : i); // Esoteric-Pull: standard store, or streaming part 1/2
		transfer_buffer[b*A+a] = fi[index]; // fpxx_copy allows direct copying without decompression+compression
	}
}
void insert_fi(const uint a, const uint A, const uint n, const uint side, const ulong t, const global fpxx_copy* transfer_buffer, global fpxx_copy* fi) {
	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	for(uint b=0u; b<def_transfers; b++) {
		const uint i = index_transfer(side*def_transfers+b);
		const ulong index = index_f(i%2u ? n : j[i-1u], t%2ul ? i : (i%2u ? i+1u : i-1u)); // Esoteric-Pull: standard load, or streaming part 2/2
		fi[index] = transfer_buffer[b*A+a]; // fpxx_copy allows direct copying without decompression+compression
	}
}
kernel void transfer_extract_fi(const uint direction, const ulong t, global uchar* transfer_buffer_p, global uchar* transfer_buffer_m, const global fpxx_copy* fi) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	extract_fi(a, A, index_extract_p(a, direction), 2u*direction+0u, t, (global fpxx_copy*) transfer_buffer_p, fi);
	extract_fi(a, A, index_extract_m(a, direction), 2u*direction+1u, t, (global fpxx_copy*) transfer_buffer_m, fi);
}
kernel void transfer__insert_fi(const uint direction, const ulong t, const global uchar* transfer_buffer_p, const global uchar* transfer_buffer_m, global fpxx_copy* fi) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	insert_fi(a, A, index_insert_p(a, direction), 2u*direction+0u, t, (const global fpxx_copy*) transfer_buffer_p, fi);
	insert_fi(a, A, index_insert_m(a, direction), 2u*direction+1u, t, (const global fpxx_copy*) transfer_buffer_m, fi);
}
// Rho, u and flags (needed if graphics are active)
void extract_rho_u_flags(const uint a, const uint A, const uint n, global char* transfer_buffer, const global float* rho, const global float* u, const global uchar* flags) {
	((global float*)transfer_buffer)[      a] = rho[               n];
	((global float*)transfer_buffer)[    A+a] = u[                 n];
	((global float*)transfer_buffer)[ 2u*A+a] = u[    def_N+(ulong)n];
	((global float*)transfer_buffer)[ 3u*A+a] = u[2ul*def_N+(ulong)n];
	((global uchar*)transfer_buffer)[16u*A+a] = flags[             n];
}
void insert_rho_u_flags(const uint a, const uint A, const uint n, const global char* transfer_buffer, global float* rho, global float* u, global uchar* flags) {
	rho[               n] = ((const global float*)transfer_buffer)[      a];
	u[                 n] = ((const global float*)transfer_buffer)[    A+a];
	u[    def_N+(ulong)n] = ((const global float*)transfer_buffer)[ 2u*A+a];
	u[2ul*def_N+(ulong)n] = ((const global float*)transfer_buffer)[ 3u*A+a];
	flags[             n] = ((const global uchar*)transfer_buffer)[16u*A+a];
}
kernel void transfer_extract_rho_u_flags(const uint direction, const ulong t, global uchar* transfer_buffer_p, global uchar* transfer_buffer_m, const global float* rho, const global float* u, const global uchar* flags) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	extract_rho_u_flags(a, A, index_extract_p(a, direction), (global char*) transfer_buffer_p, rho, u, flags);
	extract_rho_u_flags(a, A, index_extract_m(a, direction), (global char*) transfer_buffer_m, rho, u, flags);
}
kernel void transfer__insert_rho_u_flags(const uint direction, const ulong t, const global uchar* transfer_buffer_p, const global uchar* transfer_buffer_m, global float* rho, global float* u, global uchar* flags) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	insert_rho_u_flags(a, A, index_insert_p(a, direction), (const global char*) transfer_buffer_p, rho, u, flags);
	insert_rho_u_flags(a, A, index_insert_m(a, direction), (const global char*) transfer_buffer_m, rho, u, flags);
}
// Qi (Charge ddfs) 
#ifdef MAGNETO_HYDRO
void extract_gi(const uint a, const uint n, const uint side, const ulong t, global fpxx_copy* transfer_buffer, const global fpxx_copy* qi) {
	uint j7[7u]; // neighbor indices
	neighbors_charge(n, j7); // calculate neighbor indices
	const uint i = side+1u;
	const ulong index = index_f(i%2u ? j7[i] : n, t%2ul ? (i%2u ? i+1u : i-1u) : i); // Esoteric-Pull: standard store, or streaming part 1/2
	transfer_buffer[a] = qi[index]; // fpxx_copy allows direct copying without decompression+compression
}
void insert_qi(const uint a, const uint n, const uint side, const ulong t, const global fpxx_copy* transfer_buffer, global fpxx_copy* qi) {
	uint j7[7u]; // neighbor indices
	neighbors_charge(n, j7); // calculate neighbor indices
	const uint i = side+1u;
	const ulong index = index_f(i%2u ? n : j7[i-1u], t%2ul ? i : (i%2u ? i+1u : i-1u)); // Esoteric-Pull: standard load, or streaming part 2/2
	qi[index] = transfer_buffer[a]; // fpxx_copy allows direct copying without decompression+compression
}
kernel void transfer_extract_qi(const uint direction, const ulong t, global fpxx_copy* transfer_buffer_p, global fpxx_copy* transfer_buffer_m, const global fpxx_copy* qi) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	extract_gi(a, index_extract_p(a, direction), 2u*direction+0u, t, transfer_buffer_p, qi);
	extract_gi(a, index_extract_m(a, direction), 2u*direction+1u, t, transfer_buffer_m, qi);
}
kernel void transfer__insert_qi(const uint direction, const ulong t, const global fpxx_copy* transfer_buffer_p, const global fpxx_copy* transfer_buffer_m, global fpxx_copy* qi) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	insert_qi(a, index_insert_p(a, direction), 2u*direction+0u, t, transfer_buffer_p, qi);
	insert_qi(a, index_insert_m(a, direction), 2u*direction+1u, t, transfer_buffer_m, qi);
}
#endif // MAGNETO_HYDRO