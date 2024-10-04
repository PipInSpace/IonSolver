//THESE ARE TO BE REMOVED
#if defined(FP16S) || defined(FP16C)
#define fpxx ushort
#else // FP32
#define fpxx float
#endif // FP32

#define DEF_NX 20u
#define DEF_NY 20u
#define DEF_NZ 20u
#define DEF_N 8000ul

#define DEF_DX 1u
#define DEF_DY 1u
#define DEF_DZ 1u
#define DEF_DI 0u

#define DEF_AX 1u
#define DEF_AY 1u
#define DEF_AZ 1u

#define D "D2Q9" // D2Q9/D3Q15/D3Q19/D3Q27
#define DEF_VELOCITY_SET 9u // LBM velocity set (D2Q9/D3Q15/D3Q19/D3Q27)
#define DEF_DIMENSIONS 2u // number spatial dimensions (2D or 3D)
#define DEF_TRANSFERS 3u // number of DDFs that are transferred between multiple domains

#define DEF_C 0.57735027f // lattice speed of sound c = 1/sqrt(3)*dt
#define DEF_W 2.0f // relaxation rate w = dt/tau = dt/(nu/c^2+dt/2) = 1/(3*nu+1/2)
#define DEF_W0 (1.0f/2.25f)
#define DEF_WS (1.0f/9.0f)
#define DEF_WE (1.0f/36.0f)
#define DEF_KE 8.9875517923E9f
#define DEF_KMU 0.0f
#define def_ind_r 5 // Range of induction fill around cell
#define DEF_WQ 0.1f

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
#define DEF_LOD_DEPTH 2u
#define DEF_NUM_LOD 73u
#define DEF_NUM_LOD_OWN 73u
//These defines are for code completion only and are removed from the code before compilation 
#define EndTempDefines%

// Helper functions 
inline float sq(const float x) {
	return x*x;
}
inline uint uint_sq(const uint x) {
	return x*x;
}
inline float cb(const float x) {
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
float cbmagnitude(float3 v){
	return cb(sqrt(sq(v.x) + sq(v.y) + sq(v.z)));
}
inline int imax(int x, int y) {
	if (x > y) {
		return x;
	} else {
		return y;
	}
}
inline int imin(int x, int y) {
	if (x < y) {
		return x;
	} else {
		return y;
	}
}
// x to the power of DEF_DIMENSIONS
inline int to_d(int x) {
	#if defined(D2Q9)
	return x * x;
	#else
	return x * x * x;
	#endif
}
// Atomic float addition implementations for various platforms 
inline void atomic_add_f(volatile __global float* addr, const float val) {
	#if defined(cl_nv_pragma_unroll) // use hardware-supported atomic addition on Nvidia GPUs with inline PTX assembly
		float ret; asm volatile("atom.global.add.f32 %0,[%1],%2;":"=f"(ret):"l"(addr),"f"(val):"memory");
	#elif defined(__opencl_c_ext_fp32_global_atomic_add) // use hardware-supported atomic addition on some Intel GPUs
		atomic_fetch_add((volatile global atomic_float*)addr, val);
	#elif __has_builtin(__builtin_amdgcn_global_atomic_fadd_f32) // use hardware-supported atomic addition on some AMD GPUs
		__builtin_amdgcn_global_atomic_fadd_f32(addr, val);
	#else // fallback emulation: https://forums.developer.nvidia.com/t/atomicadd-float-float-atomicmul-float-float/14639/5
		float old = val; while((old=atomic_xchg(addr, atomic_xchg(addr, 0.0f)+old))!=0.0f);
	#endif
}


float3 position(const uint3 xyz) { // 3D coordinates to 3D position
	return (float3)((float)xyz.x+0.5f-0.5f*(float)DEF_NX, (float)xyz.y+0.5f-0.5f*(float)DEF_NY, (float)xyz.z+0.5f-0.5f*(float)DEF_NZ);
}
uint3 coordinates(const uint n) { // disassemble 1D index to 3D coordinates (n -> x,y,z)
	const uint t = n%(DEF_NX*DEF_NY);
	return (uint3)(t%DEF_NX, t/DEF_NX, n/(DEF_NX*DEF_NY)); // n = x+(y+z*Ny)*Nx
}
uint3 coordinates_sl(const uint n, const uint nx, const uint ny) { // disassemble 1D index and side lenghts to 3D coordinates (n -> x,y,z)
	const uint t = n%(nx*ny);
	return (uint3)(t%nx, t/nx, n/(nx*ny)); // n = x+(y+z*Ny)*Nx
}
uint index(const uint3 xyz) { // assemble 1D index from 3D coordinates (x,y,z -> n)
	return xyz.x+(xyz.y+xyz.z*DEF_NY)*DEF_NX; // n = x+(y+z*Ny)*Nx
}
bool is_halo(const uint n) {
	const uint3 xyz = coordinates(n);
	return ((DEF_DX>1u)&(xyz.x==0u||xyz.x>=DEF_NX-1u))||((DEF_DY>1u)&(xyz.y==0u||xyz.y>=DEF_NY-1u))||((DEF_DZ>1u)&(xyz.z==0u||xyz.z>=DEF_NZ-1u));
}
bool is_halo_q(const uint3 xyz) {
	return ((DEF_DX>1u)&(xyz.x==0u||xyz.x>=DEF_NX-2u))||((DEF_DY>1u)&(xyz.y==0u||xyz.y>=DEF_NY-2u))||((DEF_DZ>1u)&(xyz.z==0u||xyz.z>=DEF_NZ-2u)); // halo data is kept up-to-date, so allow using halo data for rendering
}
ulong index_f(const uint n, const uint i) { // 64-bit indexing (maximum 2^32 lattice points (1624^3 lattice resolution, 225GB)
	return (ulong)i*DEF_N+(ulong)n; // SoA (229% faster on GPU)
}
void calculate_f_eq(const float rho, float ux, float uy, float uz, float* feq) {
    const float c3=-3.0f*(sq(ux)+sq(uy)+sq(uz)), rhom1=rho-1.0f; // c3 = -2*sq(u)/(2*sq(c)), rhom1 is arithmetic optimization to minimize digit extinction
    ux *= 3.0f;
    uy *= 3.0f;
    uz *= 3.0f;
    feq[ 0] = DEF_W0*fma(rho, 0.5f*c3, rhom1); // 000 (identical for all velocity sets)
    #if defined(D2Q9)
    const float u0=ux+uy, u1=ux-uy; // these pre-calculations make manual unrolling require less FLOPs
    const float rhos=DEF_WS*rho, rhoe=DEF_WE*rho, rhom1s=DEF_WS*rhom1, rhom1e=DEF_WE*rhom1;
    feq[ 1] = fma(rhos, fma(0.5f, fma(ux, ux, c3), ux), rhom1s); feq[ 2] = fma(rhos, fma(0.5f, fma(ux, ux, c3), -ux), rhom1s); // +00 -00
    feq[ 3] = fma(rhos, fma(0.5f, fma(uy, uy, c3), uy), rhom1s); feq[ 4] = fma(rhos, fma(0.5f, fma(uy, uy, c3), -uy), rhom1s); // 0+0 0-0
    feq[ 5] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), u0), rhom1e); feq[ 6] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), -u0), rhom1e); // ++0 --0
    feq[ 7] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), u1), rhom1e); feq[ 8] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), -u1), rhom1e); // +-0 -+0
    #elif defined(D3Q15)
    const float u0=ux+uy+uz, u1=ux+uy-uz, u2=ux-uy+uz, u3=-ux+uy+uz;
    const float rhos=DEF_WS*rho, rhoc=DEF_WC*rho, rhom1s=DEF_WS*rhom1, rhom1c=DEF_WC*rhom1;
    feq[ 1] = fma(rhos, fma(0.5f, fma(ux, ux, c3), ux), rhom1s); feq[ 2] = fma(rhos, fma(0.5f, fma(ux, ux, c3), -ux), rhom1s); // +00 -00
    feq[ 3] = fma(rhos, fma(0.5f, fma(uy, uy, c3), uy), rhom1s); feq[ 4] = fma(rhos, fma(0.5f, fma(uy, uy, c3), -uy), rhom1s); // 0+0 0-0
    feq[ 5] = fma(rhos, fma(0.5f, fma(uz, uz, c3), uz), rhom1s); feq[ 6] = fma(rhos, fma(0.5f, fma(uz, uz, c3), -uz), rhom1s); // 00+ 00-
    feq[ 7] = fma(rhoc, fma(0.5f, fma(u0, u0, c3), u0), rhom1c); feq[ 8] = fma(rhoc, fma(0.5f, fma(u0, u0, c3), -u0), rhom1c); // +++ ---
    feq[ 9] = fma(rhoc, fma(0.5f, fma(u1, u1, c3), u1), rhom1c); feq[10] = fma(rhoc, fma(0.5f, fma(u1, u1, c3), -u1), rhom1c); // ++- --+
    feq[11] = fma(rhoc, fma(0.5f, fma(u2, u2, c3), u2), rhom1c); feq[12] = fma(rhoc, fma(0.5f, fma(u2, u2, c3), -u2), rhom1c); // +-+ -+-
    feq[13] = fma(rhoc, fma(0.5f, fma(u3, u3, c3), u3), rhom1c); feq[14] = fma(rhoc, fma(0.5f, fma(u3, u3, c3), -u3), rhom1c); // -++ +--
    #elif defined(D3Q19)
    const float u0=ux+uy, u1=ux+uz, u2=uy+uz, u3=ux-uy, u4=ux-uz, u5=uy-uz;
    const float rhos=DEF_WS*rho, rhoe=DEF_WE*rho, rhom1s=DEF_WS*rhom1, rhom1e=DEF_WE*rhom1;
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
    const float rhos=DEF_WS*rho, rhoe=DEF_WE*rho, rhoc=DEF_WC*rho, rhom1s=DEF_WS*rhom1, rhom1e=DEF_WE*rhom1, rhom1c=DEF_WC*rhom1;
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
    for(uint i=1u; i<DEF_VELOCITY_SET; i++) rho += f[i]; // calculate density from fi
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
	for(uint i=1u; i<DEF_VELOCITY_SET; i+=2u) {
		fhn[i   ] = load(fi, index_f(n   , t%2ul ? i    : i+1u));
		fhn[i+1u] = load(fi, index_f(j[i], t%2ul ? i+1u : i   ));
	}
}
void store_f(const uint n, const float* fhn, global fpxx* fi, const uint* j, const ulong t) {
	store(fi, index_f(n, 0u), fhn[0]); // Esoteric-Pull
	for(uint i=1u; i<DEF_VELOCITY_SET; i+=2u) {
		store(fi, index_f(j[i], t%2ul ? i+1u : i   ), fhn[i   ]);
		store(fi, index_f(n   , t%2ul ? i    : i+1u), fhn[i+1u]);
    }
}
void calculate_indices(const uint n, uint* x0, uint* xp, uint* xm, uint* y0, uint* yp, uint* ym, uint* z0, uint* zp, uint* zm) {
    const uint3 xyz = coordinates(n);
    *x0 =   xyz.x; // pre-calculate indices (periodic boundary conditions)
    *xp =  (xyz.x       +1u)%DEF_NX;
    *xm =  (xyz.x+DEF_NX-1u)%DEF_NX;
    *y0 =   xyz.y                   *DEF_NX;
    *yp = ((xyz.y       +1u)%DEF_NY)*DEF_NX;
    *ym = ((xyz.y+DEF_NY-1u)%DEF_NY)*DEF_NX;
    *z0 =   xyz.z                   *DEF_NY*DEF_NX;
    *zp = ((xyz.z       +1u)%DEF_NZ)*DEF_NY*DEF_NX;
    *zm = ((xyz.z+DEF_NZ-1u)%DEF_NZ)*DEF_NY*DEF_NX;
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
	return (float3)(u[n], u[DEF_N+(ulong)n], u[2ul*DEF_N+(ulong)n]);
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
	const float c[3u*DEF_VELOCITY_SET] = {
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
	const float w[DEF_VELOCITY_SET] = { DEF_W0, // velocity set weights
	#if defined(D2Q9)
		DEF_WS, DEF_WS, DEF_WS, DEF_WS, DEF_WE, DEF_WE, DEF_WE, DEF_WE
	#elif defined(D3Q15)
		DEF_WS, DEF_WS, DEF_WS, DEF_WS, DEF_WS, DEF_WS,
		DEF_WC, DEF_WC, DEF_WC, DEF_WC, DEF_WC, DEF_WC, DEF_WC, DEF_WC
	#elif defined(D3Q19)
		DEF_WS, DEF_WS, DEF_WS, DEF_WS, DEF_WS, DEF_WS,
		DEF_WE, DEF_WE, DEF_WE, DEF_WE, DEF_WE, DEF_WE, DEF_WE, DEF_WE, DEF_WE, DEF_WE, DEF_WE, DEF_WE
	#elif defined(D3Q27)
		DEF_WS, DEF_WS, DEF_WS, DEF_WS, DEF_WS, DEF_WS,
		DEF_WE, DEF_WE, DEF_WE, DEF_WE, DEF_WE, DEF_WE, DEF_WE, DEF_WE, DEF_WE, DEF_WE, DEF_WE, DEF_WE,
		DEF_WC, DEF_WC, DEF_WC, DEF_WC, DEF_WC, DEF_WC, DEF_WC, DEF_WC
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
	Fin[0] = 9.0f*DEF_W0*uF ; // 000 (identical for all velocity sets)
	for(uint i=1u; i<DEF_VELOCITY_SET; i++) { // loop is entirely unrolled by compiler, no unnecessary FLOPs are happening
		Fin[i] = 9.0f*w(i)*fma(c(i)*fx+c(DEF_VELOCITY_SET+i)*fy+c(2u*DEF_VELOCITY_SET+i)*fz, c(i)*ux+c(DEF_VELOCITY_SET+i)*uy+c(2u*DEF_VELOCITY_SET+i)*uz+0.33333334f, uF);
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
void load_q(const uint n, float* qhn, const global fpxx* fqi, const uint* j7, const ulong t) {
	qhn[0] = load(fqi, index_f(n, 0u)); // Esoteric-Pull
	for(uint i=1u; i<7u; i+=2u) {
		qhn[i   ] = load(fqi, index_f(n    , t%2ul ? i    : i+1u));
		qhn[i+1u] = load(fqi, index_f(j7[i], t%2ul ? i+1u : i   ));
	}
}
void store_q(const uint n, const float* qhn, global fpxx* fqi, const uint* j7, const ulong t) {
	store(fqi, index_f(n, 0u), qhn[0]); // Esoteric-Pull
	for(uint i=1u; i<7u; i+=2u) {
		store(fqi, index_f(j7[i], t%2ul ? i+1u : i   ), qhn[i   ]);
		store(fqi, index_f(n    , t%2ul ? i    : i+1u), qhn[i+1u]);
	}
}
// LOD handling
// Get 1D LOD index from cell index n and depth d (No offset for previous lod depths)
uint lod_index(const uint n, const uint d) {
	const uint3 c = coordinates(n);
	const uint nd = (1<<d); // Number of lods on each axis
	// TODO: Arithmetic optimization
	const uint x = c.x / (DEF_NX / nd);
	const uint y = c.y / (DEF_NY / nd);
	const uint z = c.z / (DEF_NZ / nd);
	return x + (y + z * nd) * nd;
}
// Size of LODs at the specified depth d
float lod_s(const uint d) {
	const uint nd = (1<<d); // Number of lods on each axis
	return (float)((DEF_NX / nd) * (DEF_NY / nd) * (DEF_NZ / nd));
}
// float coords of the center of an LOD from 1D LOD index n and specified depth d
float3 lod_coordinates(const uint n, const uint d) { // 
	const uint nd = (1<<d); // Number of lods on each axis
	const float dsx = (float)(DEF_NX / nd);
	const float dsy = (float)(DEF_NY / nd);
	const float dsz = (float)(DEF_NZ / nd);
	const uint t = n%(nd*nd);
	return (float3)((float)(t%nd)*dsx+(0.5f*dsx), (float)(t/nd)*dsy+(0.5f*dsy), (float)(n/(nd*nd))*dsz+(0.5f*dsz)); // n = x+(y+z*Ny)*Nx
}
#endif

// Simulation kernel functions
__kernel void stream_collide(global fpxx* fi, global float* rho, global float* u, global uchar* flags, const ulong t, const float fx, const float fy, const float fz 
#ifdef FORCE_FIELD
, const global float* F 
#endif // FORCE_FIELD
#ifdef MAGNETO_HYDRO
, const global float* E	// static electric field
, const global float* B	// static magnetic flux density
, global float* E_dyn	// dynamic electric field
, global float* B_dyn	// dynamic magnetic flux density
, global fpxx* fqi		// charge property of gas as ddfs
, global fpxx* ei		// electron gas ddfs
, global float* Q		// cell charge
, global float* QU_lod	// Level-of-detail for charge und velocity 
#endif // MAGNETO_HYDRO
) {
    const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
    if(n>=(uint)DEF_N||is_halo(n)) return; // don't execute stream_collide() on halo
    const uchar flagsn = flags[n]; // cache flags[n] for multiple readings
    const uchar flagsn_bo=flagsn&TYPE_BO, flagsn_su=flagsn&TYPE_SU; // extract boundary and surface flags
    if(flagsn_bo==TYPE_S||flagsn_su==TYPE_G) return; // if cell is solid boundary or gas, just return

    uint j[DEF_VELOCITY_SET]; // neighbor indices
    neighbors(n, j); // calculate neighbor indices

    float fhn[DEF_VELOCITY_SET]; // local DDFs
    load_f(n, fhn, fi, j, t); // perform streaming (part 2)

	ulong nxi=(ulong)n, nyi=DEF_N+(ulong)n, nzi=2ul*DEF_N+(ulong)n; // n indecies for x, y and z components

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

    float Fin[DEF_VELOCITY_SET]; // forcing terms, are used for ei too if MHD is enabled 
	float feq[DEF_VELOCITY_SET]; // equilibrium DDFs, are used for ei too if MHD is enabled 
    float w = DEF_W; // LBM relaxation rate w = dt/tau = dt/(nu/c^2+dt/2) = 1/(3*nu+1/2)

	#ifdef FORCE_FIELD
	{ // separate block to avoid variable name conflicts
		fxn += F[nxi]; // apply force field
		fyn += F[nyi];
		fzn += F[nzi];
	}
	#endif

	#ifdef VOLUME_FORCE
		 const float c_tau = fma(w, -0.5f, 1.0f);
	#endif // VOLUME_FORCE


	#ifdef MAGNETO_HYDRO
		/* -------- Cache fields -------- */
		float Bnx = B_dyn[nxi]; // Cache dynamic fields for multiple readings
		float Bny = B_dyn[nyi];
		float Bnz = B_dyn[nzi];
		float Enx = E_dyn[nxi];
		float Eny = E_dyn[nyi];
		float Enz = E_dyn[nzi];
		/* ---- Clear dynamic field ----- */
		B_dyn[nxi] = B[nxi]; E_dyn[nxi] = E[nxi]; // Clear dynamic buffers with static buffers for recomputation
		B_dyn[nyi] = B[nyi]; E_dyn[nyi] = E[nyi];
		B_dyn[nzi] = B[nzi]; E_dyn[nzi] = E[nzi];

		/* -------- Electron gas -------- */
		float ehn[DEF_VELOCITY_SET]; // local DDFs
		load_f(n, ehn, ei, j, t); // perform streaming (part 2)
		float e_rhon, e_uxn, e_uyn, e_uzn; // calculate local density and velocity for collision

		calculate_rho_u(ehn, &e_rhon, &e_uxn, &e_uyn, &e_uzn); // calculate (charge) density and velocity fields from ei

		// TODO: Ionization

		float e_fxn = -e_rhon * (Enx + e_uyn*Bnz - e_uzn*Bny); // F = charge * (E + (U cross B))
		float e_fyn = -e_rhon * (Eny + e_uzn*Bnx - e_uxn*Bnz); // charge is the content of the ddf
		float e_fzn = -e_rhon * (Enz + e_uxn*Bny - e_uyn*Bnx);
		const float e_rho2 = 0.5f/(e_rhon * DEF_KKGE); // convert charge density to mass density, apply external volume force (Guo forcing, Krueger p.233f)
		e_uxn = clamp(fma(fxn, e_rho2, e_uxn), -DEF_C, DEF_C); // limit velocity (for stability purposes)
		e_uyn = clamp(fma(fyn, e_rho2, e_uyn), -DEF_C, DEF_C); // force term: F*dt/(2*rho)
		e_uzn = clamp(fma(fzn, e_rho2, e_uzn), -DEF_C, DEF_C);
		calculate_forcing_terms(e_uxn, e_uyn, e_uzn, fxn, fyn, fzn, Fin); // calculate volume force terms Fin from velocity field (Guo forcing, Krueger p.233f)
		
		calculate_f_eq(e_rhon, e_uxn, e_uyn, e_uzn, feq); // calculate equilibrium DDFs

		// Perform collision with SRT TODO: use correct kinematic viscosity
		#ifdef VOLUME_FORCE
			for(uint i=0u; i<DEF_VELOCITY_SET; i++) Fin[i] *= c_tau;
		#endif // VOLUME_FORCE
		#ifndef EQUILIBRIUM_BOUNDARIES
			for(uint i=0u; i<DEF_VELOCITY_SET; i++) ehn[i] = fma(1.0f-w, ehn[i], fma(w, feq[i], Fin[i])); // perform collision (SRT)
		#else
			for(uint i=0u; i<DEF_VELOCITY_SET; i++) ehn[i] = flagsn_bo==TYPE_E ? feq[i] : fma(1.0f-w, ehn[i], fma(w, feq[i], Fin[i])); // perform collision (SRT)
		#endif // EQUILIBRIUM_BOUNDARIES

		store_f(n, ehn, ei, j, t); // perform streaming (part 1)


		/* ---- Gas charge advection ---- */
		// Advection of charge. Cell charge is stored in charge ddfs 'fqi' for advection with 'fi'.
		uint j7[7]; // neighbors of D3Q7 subset
		neighbors_charge(n, j7);
		float qhn[7]; // read from qA and stream to gh (D3Q7 subset, periodic boundary conditions)
		load_q(n, qhn, fqi, j7, t); // perform streaming (part 2)
		float Qn = 0.0f;
		for(uint i=0u; i<7u; i++) Qn += qhn[i]; // calculate charge from q
		Qn += 1.0f; // add 1.0f last to avoid digit extinction effects when summing up fqi (perturbation method / DDF-shifting)
		float qeq[7]; // cache f_equilibrium[n]
		calculate_q_eq(Qn, uxn, uyn, uzn, qeq); // calculate equilibrium DDFs
		Q[n] = Qn; // update charge field
		for(uint i=0u; i<7u; i++) qhn[i] = fma(1.0f-DEF_WQ, qhn[i], DEF_WQ*qeq[i]); // perform collision
		store_q(n, qhn, fqi, j7, t); // perform streaming (part 1)

		/* ------ EM force on gas ------- */
		// F = charge * (E + (U cross B))
		fxn += Qn * (Enx + uyn*Bnz - uzn*Bny);
		fyn += Qn * (Eny + uzn*Bnx - uxn*Bnz);
		fzn += Qn * (Enz + uxn*Bny - uyn*Bnx);


		/* ------ LOD construction ------ */
		#if DEF_LOD_DEPTH > 0 // Update LOD buffer
			uint off = 0;
			#if (DEF_DX>1 || DEF_DY>1 || DEF_DZ>1)  // Multiple Domains
				for (uint d = 0; d<=DEF_LOD_DEPTH; d++) // Iterate over depth levels d and add values to LOD buffer
			#else // Single Domain
				const uint d = DEF_LOD_DEPTH;
			#endif // Single Domain
			{
				const uint ind = (lod_index(n, d) + off) * 4;
				const float ils = 1/lod_s(d);
				atomic_add_f(&QU_lod[ind+0], Qn-e_rhon);
				atomic_add_f(&QU_lod[ind+1], uxn * ils);
				atomic_add_f(&QU_lod[ind+2], uyn * ils);
				atomic_add_f(&QU_lod[ind+3], uzn * ils);
				// offset to skip previous depths
				#if defined(D2Q9)
					off += (1<<d) * (1<<d); 
				#else
					off += (1<<d) * (1<<d) * (1<<d);
				#endif
			}
		#endif
	#endif// MAGNETO_HYDRO

	#ifdef VOLUME_FORCE
		const float rho2 = 0.5f/rhon; // apply external volume force (Guo forcing, Krueger p.233f)
		uxn = clamp(fma(fxn, rho2, uxn), -DEF_C, DEF_C); // limit velocity (for stability purposes)
		uyn = clamp(fma(fyn, rho2, uyn), -DEF_C, DEF_C); // force term: F*dt/(2*rho)
		uzn = clamp(fma(fzn, rho2, uzn), -DEF_C, DEF_C);
		calculate_forcing_terms(uxn, uyn, uzn, fxn, fyn, fzn, Fin); // calculate volume force terms Fin from velocity field (Guo forcing, Krueger p.233f)
	#else // VOLUME_FORCE
    	uxn = clamp(uxn, -DEF_C, DEF_C); // limit velocity (for stability purposes)
    	uyn = clamp(uyn, -DEF_C, DEF_C); // force term: F*dt/(2*rho)
    	uzn = clamp(uzn, -DEF_C, DEF_C);
    	for(uint i=0u; i<DEF_VELOCITY_SET; i++) Fin[i] = 0.0f;
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

    calculate_f_eq(rhon, uxn, uyn, uzn, feq); // calculate equilibrium DDFs

    #if defined(SRT) // SRT
		#ifdef VOLUME_FORCE
			for(uint i=0u; i<DEF_VELOCITY_SET; i++) Fin[i] *= c_tau;
		#endif // VOLUME_FORCE

        #ifndef EQUILIBRIUM_BOUNDARIES
            for(uint i=0u; i<DEF_VELOCITY_SET; i++) fhn[i] = fma(1.0f-w, fhn[i], fma(w, feq[i], Fin[i])); // perform collision (SRT)
        #else
            for(uint i=0u; i<DEF_VELOCITY_SET; i++) fhn[i] = flagsn_bo==TYPE_E ? feq[i] : fma(1.0f-w, fhn[i], fma(w, feq[i], Fin[i])); // perform collision (SRT)
        #endif // EQUILIBRIUM_BOUNDARIES

    #elif defined(TRT) // TRT
        const float wp = w; // TRT: inverse of "+" relaxation time
        const float wm = 1.0f/(0.1875f/(1.0f/w-0.5f)+0.5f); // TRT: inverse of "-" relaxation time wm = 1.0f/(0.1875f/(3.0f*nu)+0.5f), nu = (1.0f/w-0.5f)/3.0f;

		#ifdef VOLUME_FORCE
			const float c_taup=fma(wp, -0.25f, 0.5f), c_taum=fma(wm, -0.25f, 0.5f); // source: https://arxiv.org/pdf/1901.08766.pdf
			float Fib[DEF_VELOCITY_SET]; // F_bar
			Fib[0] = Fin[0];
			for(uint i=1u; i<DEF_VELOCITY_SET; i+=2u) {
				Fib[i   ] = Fin[i+1u];
				Fib[i+1u] = Fin[i   ];
			}
			for(uint i=0u; i<DEF_VELOCITY_SET; i++) Fin[i] = fma(c_taup, Fin[i]+Fib[i], c_taum*(Fin[i]-Fib[i]));
		#endif // VOLUME_FORCE

        float fhb[DEF_VELOCITY_SET]; // fhn in inverse directions
        float feb[DEF_VELOCITY_SET]; // feq in inverse directions
        fhb[0] = fhn[0];
        feb[0] = feq[0];
        for(uint i=1u; i<DEF_VELOCITY_SET; i+=2u) {
        	fhb[i   ] = fhn[i+1u];
        	fhb[i+1u] = fhn[i   ];
        	feb[i   ] = feq[i+1u];
        	feb[i+1u] = feq[i   ];
        }
        #ifndef EQUILIBRIUM_BOUNDARIES
            for(uint i=0u; i<DEF_VELOCITY_SET; i++) fhn[i] = fma(0.5f*wp, feq[i]-fhn[i]+feb[i]-fhb[i], fma(0.5f*wm, feq[i]-feb[i]-fhn[i]+fhb[i], fhn[i]+Fin[i])); // perform collision (TRT)
        #else // EQUILIBRIUM_BOUNDARIES
            for(uint i=0u; i<DEF_VELOCITY_SET; i++) fhn[i] = flagsn_bo==TYPE_E ? feq[i] : fma(0.5f*wp, feq[i]-fhn[i]+feb[i]-fhb[i], fma(0.5f*wm, feq[i]-feb[i]-fhn[i]+fhb[i], fhn[i]+Fin[i])); // perform collision (TRT)
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
, global fpxx* fqi
, global fpxx* ei		// electron gas ddfs
, global float* Q
#endif // MAGNETO_HYDRO
) {
    const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
    if(n>=(uint)DEF_N||is_halo(n)) return; // don't execute initialize() on halo
	ulong nxi=(ulong)n, nyi=DEF_N+(ulong)n, nzi=2ul*DEF_N+(ulong)n; // n indecies for x, y and z components
    uchar flagsn = flags[n];
    const uchar flagsn_bo = flagsn&TYPE_BO; // extract boundary flags
    uint j[DEF_VELOCITY_SET]; // neighbor indices
    neighbors(n, j); // calculate neighbor indices
    uchar flagsj[DEF_VELOCITY_SET]; // cache neighbor flags for multiple readings
    for(uint i=1u; i<DEF_VELOCITY_SET; i++) flagsj[i] = flags[j[i]];
    if(flagsn_bo==TYPE_S) { // cell is solid
	    bool TYPE_ONLY_S = true; // has only solid neighbors
	    for(uint i=1u; i<DEF_VELOCITY_SET; i++) TYPE_ONLY_S = TYPE_ONLY_S&&(flagsj[i]&TYPE_BO)==TYPE_S;
	    if(TYPE_ONLY_S) {
	    	u[nxi] = 0.0f; // reset velocity for solid lattice points with only boundary neighbors
	    	u[nyi] = 0.0f;
	    	u[nzi] = 0.0f;
			#ifdef MAGNETO_HYDRO
			Q[n] = 0.0f;
			#endif
	    }
        if(flagsn_bo==TYPE_S) {
	        u[nxi] = 0.0f; // reset velocity for all solid lattice points
	        u[nyi] = 0.0f;
	        u[nzi] = 0.0f;
			#ifdef MAGNETO_HYDRO
			Q[n] = 0.0f;
			#endif
        }
    }
    float fe_eq[DEF_VELOCITY_SET]; // f_equilibrium, reused for e_equilibrium
    calculate_f_eq(rho[n], u[n], u[DEF_N+(ulong)n], u[2ul*DEF_N+(ulong)n], fe_eq);
    store_f(n, fe_eq, fi, j, 1ul); // write to fi

	#ifdef MAGNETO_HYDRO
		// Initialize charge ddfs
		float qeq[7]; // q_equilibrium
		calculate_q_eq(Q[n], u[n], u[DEF_N+(ulong)n], u[2ul*DEF_N+(ulong)n], qeq);
		uint j7[7]; // neighbors of D3Q7 subset
		neighbors_charge(n, j7);
		store_q(n, qeq, fqi, j7, 1ul); // write to fqi. perform streaming (part 1)
		// Clear dyn with static field for recomputation
		B_dyn[nxi] = B[nxi];
		B_dyn[nyi] = B[nyi];
		B_dyn[nzi] = B[nzi];
		E_dyn[nxi] = E[nxi];
		E_dyn[nyi] = E[nyi];
		E_dyn[nzi] = E[nzi];
		// Initialize electron gas ddfs
		calculate_f_eq(Q[n], u[n], u[DEF_N+(ulong)n], u[2ul*DEF_N+(ulong)n], fe_eq);
    	store_f(n, fe_eq, ei, j, 1ul); // write to fi
	#endif // MAGNETO_HYDRO
} // initialize()

__kernel void update_fields(const global fpxx* fi, global float* rho, global float* u, const global uchar* flags, const ulong t, const float fx, const float fy, const float fz) {
    const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
    if(n>=(uint)DEF_N||is_halo(n)) return; // don't execute update_fields() on halo
    const uchar flagsn = flags[n];
    const uchar flagsn_bo=flagsn&TYPE_BO, flagsn_su=flagsn&TYPE_SU; // extract boundary and surface flags
    if(flagsn_bo==TYPE_S||flagsn_su==TYPE_G) return; // don't update fields for boundary or gas lattice points

    uint j[DEF_VELOCITY_SET]; // neighbor indices
    neighbors(n, j); // calculate neighbor indices
    float fhn[DEF_VELOCITY_SET]; // local DDFs
    load_f(n, fhn, fi, j, t); // perform streaming (part 2)

    float rhon, uxn, uyn, uzn; // calculate local density and velocity for collision
    calculate_rho_u(fhn, &rhon, &uxn, &uyn, &uzn); // calculate density and velocity fields from fi
    float fxn=fx, fyn=fy, fzn=fz; // force starts as constant volume force, can be modified before call of calculate_forcing_terms(...)
    {
        uxn = clamp(uxn, -DEF_C, DEF_C); // limit velocity (for stability purposes)
        uyn = clamp(uyn, -DEF_C, DEF_C); // force term: F*dt/(2*rho)
        uzn = clamp(uzn, -DEF_C, DEF_C);
    }

    rho[               n] = rhon; // update density field
    u[                 n] = uxn; // update velocity field
    u[    DEF_N+(ulong)n] = uyn;
    u[2ul*DEF_N+(ulong)n] = uzn;
} // update_fields()

#ifdef MAGNETO_HYDRO
__kernel void update_e_b_dynamic(global float* E_dyn, global float* B_dyn, const global float* Q, const global float* u, const global float* QU_lod, const global uchar* flags) {
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
    if(n>=(uint)DEF_N||is_halo(n)) return; // don't execute update_e_b_dynamic() on halo
    const uchar flagsn = flags[n]; // cache flags[n] for multiple readings
    const uchar flagsn_bo=flagsn&TYPE_BO, flagsn_su=flagsn&TYPE_SU; // extract boundary and surface flags
    if(flagsn_bo==TYPE_S||flagsn_su==TYPE_G) return; // if cell is solid boundary or gas, just return

	const uint3 coord_n = coordinates(n); // Cell coordinate
	const float3 coord_nf = convert_float3(coord_n); // Cell coordinate as float vector

	const uint nd = (1<<DEF_LOD_DEPTH); // Number of lowest level lods on each axis
	const uint dsx = imax(DEF_NX / (1<<nd), 1);
	const uint dsy = imax(DEF_NY / (1<<nd), 1);
	const uint dsz = imax(DEF_NZ / (1<<nd), 1);

	const uint x_upper = imin((coord_n.x / dsx) * dsx + dsx, DEF_DX>1?DEF_NX-1:DEF_NX); // Do not read at halo offsets
	const uint y_upper = imin((coord_n.y / dsy) * dsy + dsy, DEF_DY>1?DEF_NY-1:DEF_NY);
	const uint z_upper = imin((coord_n.z / dsz) * dsz + dsz, DEF_DZ>1?DEF_NZ-1:DEF_NZ);

	float3 e = {0.0f, 0.0f, 0.0f}, b = {0.0f, 0.0f, 0.0f};
	
	/// Close distance - consider individual cells
	for(		uint x = imax((coord_n.x / dsx) * dsx, DEF_DX>1?1:0); x < x_upper; x++) {
		for(	uint y = imax((coord_n.y / dsy) * dsy, DEF_DY>1?1:0); y < y_upper; y++) {
			for(uint z = imax((coord_n.z / dsz) * dsz, DEF_DZ>1?1:0); z < z_upper; z++) {
				// _c vars describe surronding cells 
				const uint n_c = x + (y + z * DEF_NY) * DEF_NX;
				if (n == n_c) continue;
					
				const float q_c = Q[n_c]; // charge of nearby cell
				if (q_c == 0.0f) { continue; } // cells without charge have no influence
				const float3 v_c = {u[n_c], u[(ulong)n_c+DEF_N], u[(ulong)n_c+DEF_N*2ul]}; // velocity of nearby cell

				// precalculation for both fields
				const float3 vec_r     = coord_nf - convert_float3(coordinates(n_c));
				const float3 pre_field = vec_r / cbmagnitude(vec_r);

				e += q_c * pre_field;             // E imparted by nearby cell (Coulomb)
				b += q_c * cross(v_c, pre_field); // B imparted by nearby cell (Biot-Savart)
			}
		}
	}

	/// Medium distance - consider lowest level LODs in own domain
	const uint ndi = lod_index(n, DEF_LOD_DEPTH); // Own LOD index, needs to be skipped
	// Loop over all lowest level LODs
	for (uint d = imax(DEF_NUM_LOD_OWN - to_d(1<<DEF_LOD_DEPTH), 0); d < DEF_NUM_LOD_OWN; d++) {
		if (d == ndi) continue;
		const float3 d_c = lod_coordinates(d, DEF_LOD_DEPTH);
		const float  q_c =  QU_lod[(d*4)+0]; // charge of LOD
		const float3 v_c = {QU_lod[(d*4)+1], QU_lod[(d*4)+2], QU_lod[(d*4)+3]}; // velocity of LOD

		// precalculation for both fields
		const float3 vec_r     = coord_nf - d_c;
		const float3 pre_field = vec_r / cbmagnitude(vec_r);

		e += q_c * pre_field;             // E imparted by LOD (Coulomb)
		b += q_c * cross(v_c, pre_field); // B imparted by LOD (Biot-Savart)
	}

	/// Large distance - consider LODs of varying detail in foreign domains (synchronized over communicate_qu_lods)
	const uint3 coord_d = coordinates_sl(DEF_DI, DEF_DX, DEF_DY); // Own domain coordinate
	uint offset = DEF_NUM_LOD_OWN;
	for (uint d = 0; d < DEF_DX*DEF_DY*DEF_DZ; d++) { // Loop over every other domain
		if (d == DEF_DI) continue;
		const uint3 coord_fd = coordinates_sl(d, DEF_DX, DEF_DY); // coordinate of foreign domain
		const int3 domain_diff = {(int)coord_d.x - (int)coord_fd.x, (int)coord_d.y - (int)coord_fd.y, (int)coord_d.z - (int)coord_fd.z}; // Difference in domain coordinates
		const uint dist = imax(abs(domain_diff.x), imax(abs(domain_diff.y), abs(domain_diff.z)));
		const uint depth = imax(0, DEF_LOD_DEPTH - dist); // Depth of foreign domain LOD
		const uint n_lod_fd = to_d(1<<depth); // Number of lods in foreign domain

		for (int l = 0; l<n_lod_fd; l++) { // Loop over every LOD in foreign domain
			float3 lc = lod_coordinates(l, depth); // LOD coordinate
			lc.x -= (float)(domain_diff.x * DEF_NX);
			lc.y -= (float)(domain_diff.y * DEF_NY);
			lc.z -= (float)(domain_diff.z * DEF_NZ);
			const float  q_c =  QU_lod[(offset+l)*4+0]; // charge of LOD
			const float3 v_c = {QU_lod[(offset+l)*4+1], QU_lod[(offset+l)*4+2], QU_lod[(offset+l)*4+3]}; // velocity of LOD

			const float3 vec_r     = coord_nf - lc;
			const float3 pre_field = vec_r / cbmagnitude(vec_r);

			e += q_c * pre_field;             // E imparted by LOD (Coulomb)
			b += q_c * cross(v_c, pre_field); // B imparted by LOD (Biot-Savart)
		}
		offset += n_lod_fd;
	}

	// update buffers
	E_dyn[n					] += DEF_KE * e.x;
	E_dyn[(ulong)n+DEF_N	] += DEF_KE * e.y;
	E_dyn[(ulong)n+DEF_N*2ul] += DEF_KE * e.z;

	B_dyn[n					] += DEF_KMU * b.x;
	B_dyn[(ulong)n+DEF_N	] += DEF_KMU * b.y;
	B_dyn[(ulong)n+DEF_N*2ul] += DEF_KMU * b.z;
} // update_e_b_dynamic()

__kernel void clear_qu_lod(global float* QU_lod) {
	// Clears own domain LODs for recomputation
	const uint n = get_global_id(0);
    if(n>DEF_NUM_LOD_OWN) return;
	QU_lod[(n * 4)+0] = 0.0f;
	QU_lod[(n * 4)+1] = 0.0f;
	QU_lod[(n * 4)+2] = 0.0f;
	QU_lod[(n * 4)+3] = 0.0f;
} // clear_qu_lod()
#endif

// Inter-Domain Transfer kernels
uint get_area(const uint direction) {
	const uint A[3] = { DEF_AX, DEF_AY, DEF_AZ };
	return A[direction];
}
// Return 1D index of cell to be transferred from id a and xyz direction (0, 1, 2) for pos and neg directions
uint index_extract_p(const uint a, const uint direction) {
	const uint3 coordinates[3] = { (uint3)(DEF_NX-2u, a%DEF_NY, a/DEF_NY), (uint3)(a/DEF_NZ, DEF_NY-2u, a%DEF_NZ), (uint3)(a%DEF_NX, a/DEF_NX, DEF_NZ-2u) };
	return index(coordinates[direction]);
}
uint index_extract_m(const uint a, const uint direction) {
	const uint3 coordinates[3] = { (uint3)(       1u, a%DEF_NY, a/DEF_NY), (uint3)(a/DEF_NZ,        1u, a%DEF_NZ), (uint3)(a%DEF_NX, a/DEF_NX,        1u) };
	return index(coordinates[direction]);
}
uint index_insert_p(const uint a, const uint direction) {
	const uint3 coordinates[3] = { (uint3)(DEF_NX-1u, a%DEF_NY, a/DEF_NY), (uint3)(a/DEF_NZ, DEF_NY-1u, a%DEF_NZ), (uint3)(a%DEF_NX, a/DEF_NX, DEF_NZ-1u) };
	return index(coordinates[direction]);
}
uint index_insert_m(const uint a, const uint direction) {
	const uint3 coordinates[3] = { (uint3)(       0u, a%DEF_NY, a/DEF_NY), (uint3)(a/DEF_NZ,        0u, a%DEF_NZ), (uint3)(a%DEF_NX, a/DEF_NX,        0u) };
	return index(coordinates[direction]);
}
// Returns an index for the transferred ddfs
uint index_transfer(const uint side_i) {
	const uchar index_transfer_data[2u*DEF_DIMENSIONS*DEF_TRANSFERS] = {
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
	uint j[DEF_VELOCITY_SET]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	for(uint b=0u; b<DEF_TRANSFERS; b++) {
		const uint i = index_transfer(side*DEF_TRANSFERS+b);
		const ulong index = index_f(i%2u ? j[i] : n, t%2ul ? (i%2u ? i+1u : i-1u) : i); // Esoteric-Pull: standard store, or streaming part 1/2
		transfer_buffer[b*A+a] = fi[index]; // fpxx_copy allows direct copying without decompression+compression
	}
}
void insert_fi(const uint a, const uint A, const uint n, const uint side, const ulong t, const global fpxx_copy* transfer_buffer, global fpxx_copy* fi) {
	uint j[DEF_VELOCITY_SET]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	for(uint b=0u; b<DEF_TRANSFERS; b++) {
		const uint i = index_transfer(side*DEF_TRANSFERS+b);
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
	((global float*)transfer_buffer)[ 2u*A+a] = u[    DEF_N+(ulong)n];
	((global float*)transfer_buffer)[ 3u*A+a] = u[2ul*DEF_N+(ulong)n];
	((global uchar*)transfer_buffer)[16u*A+a] = flags[             n];
}
void insert_rho_u_flags(const uint a, const uint A, const uint n, const global char* transfer_buffer, global float* rho, global float* u, global uchar* flags) {
	rho[               n] = ((const global float*)transfer_buffer)[      a];
	u[                 n] = ((const global float*)transfer_buffer)[    A+a];
	u[    DEF_N+(ulong)n] = ((const global float*)transfer_buffer)[ 2u*A+a];
	u[2ul*DEF_N+(ulong)n] = ((const global float*)transfer_buffer)[ 3u*A+a];
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
void extract_fqi(const uint a, const uint n, const uint side, const ulong t, global fpxx_copy* transfer_buffer, const global fpxx_copy* fqi) {
	uint j7[7u]; // neighbor indices
	neighbors_charge(n, j7); // calculate neighbor indices
	const uint i = side+1u;
	const ulong index = index_f(i%2u ? j7[i] : n, t%2ul ? (i%2u ? i+1u : i-1u) : i); // Esoteric-Pull: standard store, or streaming part 1/2
	transfer_buffer[a] = fqi[index]; // fpxx_copy allows direct copying without decompression+compression
}
void insert_fqi(const uint a, const uint n, const uint side, const ulong t, const global fpxx_copy* transfer_buffer, global fpxx_copy* fqi) {
	uint j7[7u]; // neighbor indices
	neighbors_charge(n, j7); // calculate neighbor indices
	const uint i = side+1u;
	const ulong index = index_f(i%2u ? n : j7[i-1u], t%2ul ? i : (i%2u ? i+1u : i-1u)); // Esoteric-Pull: standard load, or streaming part 2/2
	fqi[index] = transfer_buffer[a]; // fpxx_copy allows direct copying without decompression+compression
}
kernel void transfer_extract_fqi(const uint direction, const ulong t, global uchar* transfer_buffer_p, global uchar* transfer_buffer_m, const global fpxx_copy* fqi) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	extract_fqi(a, index_extract_p(a, direction), 2u*direction+0u, t, (global fpxx_copy*)transfer_buffer_p, fqi);
	extract_fqi(a, index_extract_m(a, direction), 2u*direction+1u, t, (global fpxx_copy*)transfer_buffer_m, fqi);
}
kernel void transfer__insert_fqi(const uint direction, const ulong t, const global uchar* transfer_buffer_p, const global uchar* transfer_buffer_m, global fpxx_copy* fqi) {
	const uint a=get_global_id(0), A=get_area(direction); // a = domain area index for each side, A = area of the domain boundary
	if(a>=A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
	insert_fqi(a, index_insert_p(a, direction), 2u*direction+0u, t, (const global fpxx_copy*)transfer_buffer_p, fqi);
	insert_fqi(a, index_insert_m(a, direction), 2u*direction+1u, t, (const global fpxx_copy*)transfer_buffer_m, fqi);
}
#endif // MAGNETO_HYDRO
