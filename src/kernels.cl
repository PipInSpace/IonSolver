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

#define def_Dx 1u
#define def_Dy 1u
#define def_Dz 1u

#define def_Ox 1 // offsets are signed integer!
#define def_Oy 1
#define def_Oz 1

#define def_Ax 1u
#define def_Ay 1u
#define def_Az 1u

#define def_domain_offset_x 0.0f
#define def_domain_offset_y 0.0f
#define def_domain_offset_z 0.0f

#define D "D2Q9" // D2Q9/D3Q15/D3Q19/D3Q27
#define def_velocity_set 9u // LBM velocity set (D2Q9/D3Q15/D3Q19/D3Q27)
#define def_dimensions 2u // number spatial dimensions (2D or 3D)
#define def_transfers 3u // number of DDFs that are transferred between multiple domains

#define def_c 0.57735027f // lattice speed of sound c = 1/sqrt(3)*dt
#define def_w 2.0f // relaxation rate w = dt/tau = dt/(nu/c^2+dt/2) = 1/(3*nu+1/2)
#define def_w0 (1.0f/2.25f)
#define def_ws (1.0f/9.0f)
#define def_we (1.0f/36.0f)

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

#define GRAPHICS
#define def_streamline_sparse 4u
#define def_streamline_length 128u
#define def_screen_width 1920u
#define def_screen_height 1080u
#define def_scale_u 1.0f
#define def_scale_Q_min 0.0001f

#define COLOR_S (127<<16|127<<8|127)
#define COLOR_E (  0<<16|255<<8|  0)
#define COLOR_M (255<<16|  0<<8|255)
#define COLOR_T (255<<16|  0<<8|  0)
#define COLOR_F (  0<<16|  0<<8|255)
#define COLOR_I (  0<<16|255<<8|255)
#define COLOR_0 (127<<16|127<<8|127)
#define COLOR_X (255<<16|127<<8|  0)
#define COLOR_Y (255<<16|255<<8|  0)
#define COLOR_P (255<<16|255<<8|191)

//These defines are for code completion only and are removed from the code before compilation 
#define EndTempDefines%

//Graphics Helper functions:
#ifdef GRAPHICS
int color_average(const int c1, const int c2) { // (c1+c2)/s
	const uchar4 cc1=as_uchar4(c1), cc2=as_uchar4(c2);
	return as_int((uchar4)((uchar)((cc1.x+cc2.x)/2u), (uchar)((cc1.y+cc2.y)/2u), (uchar)((cc1.z+cc2.z)/2u), (uchar)0u));
}
int shading(const int c, const float3 p, const float3 normal, const float* camera_cache) {
    const float dis  = camera_cache[ 1]; // fetch camera parameters (rotation matrix, camera position, etc.)
    const float posx = camera_cache[ 2]-def_domain_offset_x;
    const float posy = camera_cache[ 3]-def_domain_offset_y;
    const float posz = camera_cache[ 4]-def_domain_offset_z;
    const float Rzx  = camera_cache[11];
    const float Rzy  = camera_cache[12];
    const float Rzz  = camera_cache[13];
    const uchar cr=c>>16&255, cg=c>>8&255, cb=c&255;
    const float nl2 = sq(normal.x)+sq(normal.y)+sq(normal.z); // only one native_sqrt instead of two
    const float dx = p.x-fma(Rzx, dis, posx); // direction of light source is viewing direction
    const float dy = p.y-fma(Rzy, dis, posy);
    const float dz = p.z-fma(Rzz, dis, posz);
    const float dl2 = sq(dx)+sq(dy)+sq(dz);
    const float br = max(1.5f*fabs(normal.x*dx+normal.y*dy+normal.z*dz)*rsqrt(nl2*dl2), 0.3f);
    return min((int)(br*cr), 255)<<16|min((int)(br*cg), 255)<<8|min((int)(br*cb), 255);
}
void draw(const int x, const int y, const float z, const int color, global int* bitmap, volatile global int* zbuffer, const int stereo) {
	const int index=x+y*def_screen_width, iz=(int)(z*(2147483647.0f/10000.0f)); // use int z-buffer and atomic_max to minimize noise in image
}
void convert_line(const float3 p0, const float3 p1, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer, const int stereo) { // 3D -> 2D
	int r0x, r0y, r1x, r1y; float r0z, r1z;
	if(convert(&r0x, &r0y, &r0z, p0, camera_cache, stereo) && convert(&r1x, &r1y, &r1z, p1, camera_cache, stereo)
		&& !(is_off_screen(r0x, r0y, stereo) && is_off_screen(r1x, r1y, stereo))) { // cancel drawing if both points are off screen
		int x=r0x, y=r0y; // Bresenham algorithm
		const float z = 0.5f*(r0z+r1z); // approximate line z position for each pixel to be equal
		const int dx= abs(r1x-r0x), sx=2*(r0x<r1x)-1;
		const int dy=-abs(r1y-r0y), sy=2*(r0y<r1y)-1;
		int err = dx+dy;
		while(x!=r1x||y!=r1y) {
			draw(x, y, z, color, bitmap, zbuffer, stereo);
			const int e2 = 2*err;
			if(e2>dy) { err+=dy; x+=sx; }
			if(e2<dx) { err+=dx; y+=sy; }
		}
	}
}
void draw_line(const float3 p0, const float3 p1, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer) { // 3D -> 2D
	const bool vr = (as_int(camera_cache[14])>>31)&0x1;
	if(!vr) {
		convert_line(p0, p1, color, camera_cache, bitmap, zbuffer,  0);
	} else {
		convert_line(p0, p1, color, camera_cache, bitmap, zbuffer, -1); // left eye
		convert_line(p0, p1, color, camera_cache, bitmap, zbuffer, +1); // right eye
	}
}
uint marching_cubes(const float* v, const float iso, float3* triangles) { // input: 8 values v, isovalue; output: returns number of triangles, 15 triangle vertices t
	uint cube = 0u; // determine index of which vertices are inside of the isosurface
	for(uint i=0u; i<8u; i++) cube |= (v[i]<iso)<<i;
	if(cube==0u||cube==255u) return 0u; // cube is entirely inside/outside of the isosurface
	float3 p[8]; // definition of unit cube corners
	p[0] = (float3)(0.0f, 0.0f, 0.0f);
	p[1] = (float3)(1.0f, 0.0f, 0.0f);
	p[2] = (float3)(1.0f, 0.0f, 1.0f);
	p[3] = (float3)(0.0f, 0.0f, 1.0f);
	p[4] = (float3)(0.0f, 1.0f, 0.0f);
	p[5] = (float3)(1.0f, 1.0f, 0.0f);
	p[6] = (float3)(1.0f, 1.0f, 1.0f);
	p[7] = (float3)(0.0f, 1.0f, 1.0f);
	const uint edges = edge_table(cube);
	float3 vertex[12]; // find the vertices where the surface intersects the cube
	if(edges&   1u) vertex[ 0] = interpolate_vertex(p[0], p[1], v[0], v[1], iso); // calculate vertices on all 12 edges
	if(edges&   2u) vertex[ 1] = interpolate_vertex(p[1], p[2], v[1], v[2], iso);
	if(edges&   4u) vertex[ 2] = interpolate_vertex(p[2], p[3], v[2], v[3], iso);
	if(edges&   8u) vertex[ 3] = interpolate_vertex(p[3], p[0], v[3], v[0], iso);
	if(edges&  16u) vertex[ 4] = interpolate_vertex(p[4], p[5], v[4], v[5], iso);
	if(edges&  32u) vertex[ 5] = interpolate_vertex(p[5], p[6], v[5], v[6], iso);
	if(edges&  64u) vertex[ 6] = interpolate_vertex(p[6], p[7], v[6], v[7], iso);
	if(edges& 128u) vertex[ 7] = interpolate_vertex(p[7], p[4], v[7], v[4], iso);
	if(edges& 256u) vertex[ 8] = interpolate_vertex(p[0], p[4], v[0], v[4], iso);
	if(edges& 512u) vertex[ 9] = interpolate_vertex(p[1], p[5], v[1], v[5], iso);
	if(edges&1024u) vertex[10] = interpolate_vertex(p[2], p[6], v[2], v[6], iso);
	if(edges&2048u) vertex[11] = interpolate_vertex(p[3], p[7], v[3], v[7], iso);
	cube *= 15u;
	uint i; // number of triangle vertices
	for(i=0u; i<15u&&triangle_table(cube+i)!=15u; i+=3u) { // create the triangles
		triangles[i   ] = vertex[triangle_table(cube+i   )];
		triangles[i+1u] = vertex[triangle_table(cube+i+1u)];
		triangles[i+2u] = vertex[triangle_table(cube+i+2u)];
	}
	return i/3u; // return number of triangles
}


#endif

float3 position(const uint3 xyz) { // 3D coordinates to 3D position
	return (float3)((float)xyz.x+0.5f-0.5f*(float)def_Nx, (float)xyz.y+0.5f-0.5f*(float)def_Ny, (float)xyz.z+0.5f-0.5f*(float)def_Nz);
}

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

    #ifndef EQUILIBRIUM_BOUNDARIES // EQUILIBRIUM_BOUNDARIES
        calculate_rho_u(fhn, &rhon, &uxn, &uyn, &uzn); // calculate density and velocity fields from fi
    #else
        if(flagsn_bo==TYPE_E) {
        	rhon = rho[               n]; // apply preset velocity/density
        	uxn  = u[                 n];
        	uyn  = u[    def_N+(ulong)n];
        	uzn  = u[2ul*def_N+(ulong)n];
        } else {
        	calculate_rho_u(fhn, &rhon, &uxn, &uyn, &uzn); // calculate density and velocity fields from fi
        }
    #endif // EQUILIBRIUM_BOUNDARIES

    float fxn=fx, fyn=fy, fzn=fz; // force starts as constant volume force, can be modified before call of calculate_forcing_terms(...)
    float Fin[def_velocity_set]; // forcing terms

    uxn = clamp(uxn, -def_c, def_c); // limit velocity (for stability purposes)
    uyn = clamp(uyn, -def_c, def_c); // force term: F*dt/(2*rho)
    uzn = clamp(uzn, -def_c, def_c);
    for(uint i=0u; i<def_velocity_set; i++) Fin[i] = 0.0f;

    float feq[def_velocity_set]; // equilibrium DDFs
    calculate_f_eq(rhon, uxn, uyn, uzn, feq); // calculate equilibrium DDFs
    float w = def_w; // LBM relaxation rate w = dt/tau = dt/(nu/c^2+dt/2) = 1/(3*nu+1/2)

    #if defined(SRT) // SRT
        #ifndef EQUILIBRIUM_BOUNDARIES
            for(uint i=0u; i<def_velocity_set; i++) fhn[i] = fma(1.0f-w, fhn[i], fma(w, feq[i], Fin[i])); // perform collision (SRT)
        #else
            for(uint i=0u; i<def_velocity_set; i++) fhn[i] = flagsn_bo==TYPE_E ? feq[i] : fma(1.0f-w, fhn[i], fma(w, feq[i], Fin[i])); // perform collision (SRT)
        #endif // EQUILIBRIUM_BOUNDARIES
    #elif defined(TRT) // TRT
        const float wp = w; // TRT: inverse of "+" relaxation time
        const float wm = 1.0f/(0.1875f/(1.0f/w-0.5f)+0.5f); // TRT: inverse of "-" relaxation time wm = 1.0f/(0.1875f/(3.0f*nu)+0.5f), nu = (1.0f/w-0.5f)/3.0f;

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

// Graphics code
#ifdef GRAPHICS
kernel void graphics_flags(const global uchar* flags, const global float* camera, global int* bitmap, global int* zbuffer) {
    const uint n = get_global_id(0);
    if(n>=(uint)def_N||is_halo(n)) return; // don't execute graphics_flags() on halo
    const uchar flagsn = flags[n]; // cache flags
    const uchar flagsn_bo = flagsn&TYPE_BO; // extract boundary flags
    if(flagsn==0u||flagsn==TYPE_G) return; // don't draw regular fluid cells
    //if(flagsn&TYPE_SU) return; // don't draw surface
    float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
    for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
    uint x0, xp, xm, y0, yp, ym, z0, zp, zm;
    calculate_indices(n, &x0, &xp, &xm, &y0, &yp, &ym, &z0, &zp, &zm);
    const uint3 xyz = coordinates(n);
    const float3 p = position(xyz);
    const int c =  // coloring scheme
    	flagsn_bo==TYPE_S ? COLOR_S : // solid boundary
    	((flagsn&TYPE_T)&&flagsn_bo==TYPE_E) ? color_average(COLOR_T, COLOR_E) : // both temperature boundary and equilibrium boundary
    	((flagsn&TYPE_T)&&flagsn_bo==TYPE_MS) ? color_average(COLOR_T, COLOR_M) : // both temperature boundary and moving boundary
    	flagsn&TYPE_T ? COLOR_T : // temperature boundary
    	flagsn_bo==TYPE_E ? COLOR_E : // equilibrium boundary
    	flagsn_bo==TYPE_MS ? COLOR_M : // moving boundary
    	flagsn&TYPE_F ? COLOR_F : // fluid
    	flagsn&TYPE_I ? COLOR_I : // interface
    	flagsn&TYPE_X ? COLOR_X : // reserved type X
    	flagsn&TYPE_Y ? COLOR_Y : // reserved type Y
    	COLOR_0; // regular or gas cell
    //draw_point(p, c, camera_cache, bitmap, zbuffer); // draw one pixel for every boundary cell
    uint t;
    t = xp+y0+z0; const bool not_xp = xyz.x<def_Nx-1u && flagsn==flags[t] && !is_halo(t); // +00
    t = xm+y0+z0; const bool not_xm = xyz.x>       0u && flagsn==flags[t] && !is_halo(t); // -00
    t = x0+yp+z0; const bool not_yp = xyz.y<def_Ny-1u && flagsn==flags[t] && !is_halo(t); // 0+0
    t = x0+ym+z0; const bool not_ym = xyz.y>       0u && flagsn==flags[t] && !is_halo(t); // 0-0
    t = x0+y0+zp; const bool not_zp = xyz.z<def_Nz-1u && flagsn==flags[t] && !is_halo(t); // 00+
    t = x0+y0+zm; const bool not_zm = xyz.z>       0u && flagsn==flags[t] && !is_halo(t); // 00-
    const float3 p0 = (float3)(p.x-0.5f, p.y-0.5f, p.z-0.5f); // ---
    const float3 p1 = (float3)(p.x+0.5f, p.y+0.5f, p.z+0.5f); // +++
    const float3 p2 = (float3)(p.x-0.5f, p.y-0.5f, p.z+0.5f); // --+
    const float3 p3 = (float3)(p.x+0.5f, p.y+0.5f, p.z-0.5f); // ++-
    const float3 p4 = (float3)(p.x-0.5f, p.y+0.5f, p.z-0.5f); // -+-
    const float3 p5 = (float3)(p.x+0.5f, p.y-0.5f, p.z+0.5f); // +-+
    const float3 p6 = (float3)(p.x+0.5f, p.y-0.5f, p.z-0.5f); // +--
    const float3 p7 = (float3)(p.x-0.5f, p.y+0.5f, p.z+0.5f); // -++
    if(!(not_xm||not_ym)) draw_line(p0, p2, c, camera_cache, bitmap, zbuffer); // to draw the entire surface, replace || by &&
    if(!(not_xm||not_zm)) draw_line(p0, p4, c, camera_cache, bitmap, zbuffer);
    if(!(not_ym||not_zm)) draw_line(p0, p6, c, camera_cache, bitmap, zbuffer);
    if(!(not_xp||not_yp)) draw_line(p1, p3, c, camera_cache, bitmap, zbuffer);
    if(!(not_xp||not_zp)) draw_line(p1, p5, c, camera_cache, bitmap, zbuffer);
    if(!(not_yp||not_zp)) draw_line(p1, p7, c, camera_cache, bitmap, zbuffer);
    if(!(not_ym||not_zp)) draw_line(p2, p5, c, camera_cache, bitmap, zbuffer);
    if(!(not_xm||not_zp)) draw_line(p2, p7, c, camera_cache, bitmap, zbuffer);
    if(!(not_yp||not_zm)) draw_line(p3, p4, c, camera_cache, bitmap, zbuffer);
    if(!(not_xp||not_zm)) draw_line(p3, p6, c, camera_cache, bitmap, zbuffer);
    if(!(not_xm||not_yp)) draw_line(p4, p7, c, camera_cache, bitmap, zbuffer);
    if(!(not_xp||not_ym)) draw_line(p5, p6, c, camera_cache, bitmap, zbuffer);
}

kernel void graphics_flags_mc(const global uchar* flags, const global float* camera, global int* bitmap, global int* zbuffer) {
    const uint n = get_global_id(0);
    if(n>=(uint)def_N||is_halo(n)) return; // don't execute graphics_flags() on halo
    const uint3 xyz = coordinates(n);
    if(xyz.x>=def_Nx-1u||xyz.y>=def_Ny-1u||xyz.z>=def_Nz-1u) return;
    //if(xyz.x==0u||xyz.y==0u||xyz.z==0u||xyz.x>=def_Nx-2u||xyz.y>=def_Ny-2u||xyz.z>=def_Nz-2u) return;
    uint j[8];
    const uint x0 =  xyz.x; // cube stencil
    const uint xp =  xyz.x+1u;
    const uint y0 =  xyz.y    *def_Nx;
    const uint yp = (xyz.y+1u)*def_Nx;
    const uint z0 =  xyz.z    *def_Ny*def_Nx;
    const uint zp = (xyz.z+1u)*def_Ny*def_Nx;
    j[0] = n       ; // 000
    j[1] = xp+y0+z0; // +00
    j[2] = xp+y0+zp; // +0+
    j[3] = x0+y0+zp; // 00+
    j[4] = x0+yp+z0; // 0+0
    j[5] = xp+yp+z0; // ++0
    j[6] = xp+yp+zp; // +++
    j[7] = x0+yp+zp; // 0++
    float v[8];
    for(uint i=0u; i<8u; i++) v[i] = (float)((flags[j[i]]&TYPE_BO)==TYPE_S);
    float3 triangles[15]; // maximum of 5 triangles with 3 vertices each
    const uint tn = marching_cubes(v, 0.5f, triangles); // run marching cubes algorithm
    if(tn==0u) return;
    float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
    for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
    const float3 offset = (float3)((float)xyz.x+0.5f-0.5f*(float)def_Nx, (float)xyz.y+0.5f-0.5f*(float)def_Ny, (float)xyz.z+0.5f-0.5f*(float)def_Nz);
    for(uint i=0u; i<tn; i++) {
	    const float3 p0 = triangles[3u*i   ];
	    const float3 p1 = triangles[3u*i+1u];
	    const float3 p2 = triangles[3u*i+2u];
	    const float3 normal = normalize(cross(p1-p0, p2-p0));

        const int c = shading(0xDFDFDF, (p0+p1+p2)/3.0f+offset, normal, camera_cache); // 0xDFDFDF;
        draw_triangle(p0+offset, p1+offset, p2+offset, c, camera_cache, bitmap, zbuffer);
    }
}

kernel void graphics_field(const global uchar* flags, const global float* u, const global float* camera, global int* bitmap, global int* zbuffer, const int slice_mode, const int slice_x, const int slice_y, const int slice_z) {
	const uint n = get_global_id(0);
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute graphics_field() on halo
	const uint3 xyz = coordinates(n);
	const bool rx=(int)xyz.x!=slice_x, ry=(int)xyz.y!=slice_y, rz=(int)xyz.z!=slice_z;
	if((slice_mode==1&&rx)||(slice_mode==2&&ry)||(slice_mode==3&&rz)||(slice_mode==4&&rx&&rz)||(slice_mode==5&&rx&&ry&&rz)||(slice_mode==6&&ry&&rz)||(slice_mode==7&&rx&&ry)) return;

    if(flags[n]&(TYPE_S|TYPE_E|TYPE_I|TYPE_G)) return;

    float3 un = load_u(n, u); // cache velocity
    const float ul = length(un);
    if(def_scale_u*ul<0.1f) return; // don't draw lattice points where the velocity is lower than this threshold
    float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
    for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
    const float3 p = position(coordinates(n));
    const int c = iron_color(255.0f*def_scale_u*ul); // coloring by velocity
    draw_line(p-(0.5f/ul)*un, p+(0.5f/ul)*un, c, camera_cache, bitmap, zbuffer);

}

kernel void graphics_streamline(const global uchar* flags, const global float* u, const global float* camera, global int* bitmap, global int* zbuffer, const int slice_mode, const int slice_x, const int slice_y, const int slice_z) {
    const uint n = get_global_id(0);
	const float3 slice = position((uint3)(slice_x, slice_y, slice_z));
    #ifndef D2Q9
        if(n>=(def_Nx/def_streamline_sparse)*(def_Ny/def_streamline_sparse)*(def_Nz/def_streamline_sparse)) return;
        const uint z = n/((def_Nx/def_streamline_sparse)*(def_Ny/def_streamline_sparse)); // disassemble 1D index to 3D coordinates
        const uint t = n%((def_Nx/def_streamline_sparse)*(def_Ny/def_streamline_sparse));
        const uint y = t/(def_Nx/def_streamline_sparse);
        const uint x = t%(def_Nx/def_streamline_sparse);
        float3 p = (float)def_streamline_sparse*((float3)((float)x+0.5f, (float)y+0.5f, (float)z+0.5f))-0.5f*((float3)((float)def_Nx, (float)def_Ny, (float)def_Nz));
        const bool rx=fabs(p.x-slice.x)>0.5f*(float)def_streamline_sparse, ry=fabs(p.y-slice.y)>0.5f*(float)def_streamline_sparse, rz=fabs(p.z-slice.z)>0.5f*(float)def_streamline_sparse;
    #else
        if(n>=(def_Nx/def_streamline_sparse)*(def_Ny/def_streamline_sparse)) return;
        const uint y = n/(def_Nx/def_streamline_sparse); // disassemble 1D index to 3D coordinates
        const uint x = n%(def_Nx/def_streamline_sparse);
        float3 p = ((float3)((float)def_streamline_sparse*((float)x+0.5f), (float)def_streamline_sparse*((float)y+0.5f), 0.5f))-0.5f*((float3)((float)def_Nx, (float)def_Ny, (float)def_Nz));
        const bool rx=fabs(p.x-slice.x)>0.5f*(float)def_streamline_sparse, ry=fabs(p.y-slice.y)>0.5f*(float)def_streamline_sparse, rz=true;
    #endif
    if((slice_mode==1&&rx)||(slice_mode==2&&ry)||(slice_mode==3&&rz)||(slice_mode==4&&rx&&rz)||(slice_mode==5&&rx&&ry&&rz)||(slice_mode==6&&ry&&rz)||(slice_mode==7&&rx&&ry)) return;
    if((slice_mode==1||slice_mode==5||slice_mode==4||slice_mode==7)&!rx) p.x = slice.x; // snap streamline position to slice position
    if((slice_mode==2||slice_mode==5||slice_mode==6||slice_mode==7)&!ry) p.y = slice.y;
    if((slice_mode==3||slice_mode==5||slice_mode==4||slice_mode==6)&!rz) p.z = slice.z;
    float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
    for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
    const float hLx=0.5f*(float)(def_Nx-2u*(def_Dx>1u)), hLy=0.5f*(float)(def_Ny-2u*(def_Dy>1u)), hLz=0.5f*(float)(def_Nz-2u*(def_Dz>1u));
    //draw_circle(p, 0.5f*def_streamline_sparse, 0xFFFFFF, camera_cache, bitmap, zbuffer);
    for(float dt=-1.0f; dt<=1.0f; dt+=2.0f) { // integrate forward and backward in time
    	float3 p0, p1=p;
    	for(uint l=0u; l<def_streamline_length/2u; l++) {
    		const uint x = (uint)(p1.x+1.5f*(float)def_Nx)%def_Nx;
    		const uint y = (uint)(p1.y+1.5f*(float)def_Ny)%def_Ny;
    		const uint z = (uint)(p1.z+1.5f*(float)def_Nz)%def_Nz;
    		const uint n = x+(y+z*def_Ny)*def_Nx;
    		if(flags[n]&(TYPE_S|TYPE_E|TYPE_I|TYPE_G)) return;
    		const float3 un = load_u(n, u); // interpolate_u(p1, u)
    		const float ul = length(un);
    		p0 = p1;
    		p1 += (dt/ul)*un; // integrate forward in time
    		if(def_scale_u*ul<0.1f||p1.x<-hLx||p1.x>hLx||p1.y<-hLy||p1.y>hLy||p1.z<-hLz||p1.z>hLz) break;
            const int c = iron_color(255.0f*def_scale_u*ul);
            draw_line(p0, p1, c, camera_cache, bitmap, zbuffer);
        }
    }
}

kernel void graphics_q_field(const global uchar* flags, const global float* u, const global float* camera, global int* bitmap, global int* zbuffer) {
    const uint n = get_global_id(0);
    if(n>=(uint)def_N||is_halo(n)) return; // don't execute graphics_q_field() on halo
    if(flags[n]&(TYPE_S|TYPE_E|TYPE_I|TYPE_G)) return;
    float3 un = load_u(n, u); // cache velocity
    const float ul = length(un);
    const float Q = calculate_Q(n, u);
    if(Q<def_scale_Q_min||ul==0.0f) return; // don't draw lattice points where the velocity is very low
    float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
    for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
    const float3 p = position(coordinates(n));
    const int c = rainbow_color(255.0f*def_scale_u*ul); // coloring by velocity
    draw_line(p-(0.5f/ul)*un, p+(0.5f/ul)*un, c, camera_cache, bitmap, zbuffer);
}

kernel void graphics_q(const global uchar* flags, const global float* u, const global float* camera, global int* bitmap, global int* zbuffer) {
    const uint n = get_global_id(0);
	const uint3 xyz = coordinates(n);
	if(xyz.x>=def_Nx-1u||xyz.y>=def_Ny-1u||xyz.z>=def_Nz-1u||is_halo_q(xyz)) return; // don't execute graphics_q_field() on marching-cubes halo
	const uint x0 =  xyz.x; // cube stencil
	const uint xp =  xyz.x+1u;
	const uint y0 =  xyz.y    *def_Nx;
	const uint yp = (xyz.y+1u)*def_Nx;
	const uint z0 =  xyz.z    *def_Ny*def_Nx;
	const uint zp = (xyz.z+1u)*def_Ny*def_Nx;
	const uint xq =  (xyz.x       +2u)%def_Nx; // central difference stencil on each cube corner point
	const uint xm =  (xyz.x+def_Nx-1u)%def_Nx;
	const uint yq = ((xyz.y       +2u)%def_Ny)*def_Nx;
	const uint ym = ((xyz.y+def_Ny-1u)%def_Ny)*def_Nx;
	const uint zq = ((xyz.z       +2u)%def_Nz)*def_Ny*def_Nx;
	const uint zm = ((xyz.z+def_Nz-1u)%def_Nz)*def_Ny*def_Nx;
	uint j[32];
	j[ 0] = n       ; // 000 // cube stencil
	j[ 1] = xp+y0+z0; // +00
	j[ 2] = xp+y0+zp; // +0+
	j[ 3] = x0+y0+zp; // 00+
	j[ 4] = x0+yp+z0; // 0+0
	j[ 5] = xp+yp+z0; // ++0
	j[ 6] = xp+yp+zp; // +++
	j[ 7] = x0+yp+zp; // 0++
	j[ 8] = xm+y0+z0; // -00 // central difference stencil on each cube corner point
	j[ 9] = x0+ym+z0; // 0-0
	j[10] = x0+y0+zm; // 00-
	j[11] = xq+y0+z0; // #00
	j[12] = xp+ym+z0; // +-0
	j[13] = xp+y0+zm; // +0-
	j[14] = xq+y0+zp; // #0+
	j[15] = xp+ym+zp; // +-+
	j[16] = xp+y0+zq; // +0#
	j[17] = xm+y0+zp; // -0+
	j[18] = x0+ym+zp; // 0-+
	j[19] = x0+y0+zq; // 00#
	j[20] = xm+yp+z0; // -+0
	j[21] = x0+yq+z0; // 0#0
	j[22] = x0+yp+zm; // 0+-
	j[23] = xq+yp+z0; // #+0
	j[24] = xp+yq+z0; // +#0
	j[25] = xp+yp+zm; // ++-
	j[26] = xq+yp+zp; // #++
	j[27] = xp+yq+zp; // +#+
	j[28] = xp+yp+zq; // ++#
	j[29] = xm+yp+zp; // -++
	j[30] = x0+yq+zp; // 0#+
	j[31] = x0+yp+zq; // 0+#
	uchar flags_cell = 0u;
	for(uint i=0u; i<32u; i++) flags_cell |= flags[j[i]];
	if(flags_cell&(TYPE_S|TYPE_E|TYPE_I|TYPE_G)) return;
	float3 uj[8];
	for(uint i=0u; i<8u; i++) uj[i] = load_u(j[i], u);
	float v[8]; // don't load any velocity twice from global memory
	v[0] = calculate_Q_cached(      uj[ 1]    , load_u(j[ 8], u),       uj[ 4]    , load_u(j[ 9], u),       uj[ 3]    , load_u(j[10], u));
	v[1] = calculate_Q_cached(load_u(j[11], u),       uj[ 0]    ,       uj[ 5]    , load_u(j[12], u),       uj[ 2]    , load_u(j[13], u));
	v[2] = calculate_Q_cached(load_u(j[14], u),       uj[ 3]    ,       uj[ 6]    , load_u(j[15], u), load_u(j[16], u),       uj[ 1]    );
	v[3] = calculate_Q_cached(      uj[ 2]    , load_u(j[17], u),       uj[ 7]    , load_u(j[18], u), load_u(j[19], u),       uj[ 0]    );
	v[4] = calculate_Q_cached(      uj[ 5]    , load_u(j[20], u), load_u(j[21], u),       uj[ 0]    ,       uj[ 7]    , load_u(j[22], u));
	v[5] = calculate_Q_cached(load_u(j[23], u),       uj[ 4]    , load_u(j[24], u),       uj[ 1]    ,       uj[ 6]    , load_u(j[25], u));
	v[6] = calculate_Q_cached(load_u(j[26], u),       uj[ 7]    , load_u(j[27], u),       uj[ 2]    , load_u(j[28], u),       uj[ 5]    );
	v[7] = calculate_Q_cached(      uj[ 6]    , load_u(j[29], u), load_u(j[30], u),       uj[ 3]    , load_u(j[31], u),       uj[ 4]    );
	float3 triangles[15]; // maximum of 5 triangles with 3 vertices each
	const uint tn = marching_cubes(v, def_scale_Q_min, triangles); // run marching cubes algorithm
	if(tn==0u) return;
	float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
	for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
	const float3 offset = (float3)((float)xyz.x+0.5f-0.5f*(float)def_Nx, (float)xyz.y+0.5f-0.5f*(float)def_Ny, (float)xyz.z+0.5f-0.5f*(float)def_Nz);
	for(uint i=0u; i<tn; i++) {
		const float3 p0 = triangles[3u*i   ]; // triangle coordinates in [0,1] (local cell)
		const float3 p1 = triangles[3u*i+1u];
		const float3 p2 = triangles[3u*i+2u];
		const float3 normal = normalize(cross(p1-p0, p2-p0));
		int c0, c1, c2; {
			const float x1=p0.x, y1=p0.y, z1=p0.z, x0=1.0f-x1, y0=1.0f-y1, z0=1.0f-z1; // calculate interpolation factors
			const float3 ui = (x0*y0*z0)*uj[0]+(x1*y0*z0)*uj[1]+(x1*y0*z1)*uj[2]+(x0*y0*z1)*uj[3]+(x0*y1*z0)*uj[4]+(x1*y1*z0)*uj[5]+(x1*y1*z1)*uj[6]+(x0*y1*z1)*uj[7]; // perform trilinear interpolation
			c0 = shading(rainbow_color(255.0f*def_scale_u*length(ui)), p0+offset, normal, camera_cache); // rainbow_color(255.0f*def_scale_u*length(ui));
		} {
			const float x1=p1.x, y1=p1.y, z1=p1.z, x0=1.0f-x1, y0=1.0f-y1, z0=1.0f-z1; // calculate interpolation factors
			const float3 ui = (x0*y0*z0)*uj[0]+(x1*y0*z0)*uj[1]+(x1*y0*z1)*uj[2]+(x0*y0*z1)*uj[3]+(x0*y1*z0)*uj[4]+(x1*y1*z0)*uj[5]+(x1*y1*z1)*uj[6]+(x0*y1*z1)*uj[7]; // perform trilinear interpolation
			c1 = shading(rainbow_color(255.0f*def_scale_u*length(ui)), p1+offset, normal, camera_cache); // rainbow_color(255.0f*def_scale_u*length(ui));
		} {
			const float x1=p2.x, y1=p2.y, z1=p2.z, x0=1.0f-x1, y0=1.0f-y1, z0=1.0f-z1; // calculate interpolation factors
			const float3 ui = (x0*y0*z0)*uj[0]+(x1*y0*z0)*uj[1]+(x1*y0*z1)*uj[2]+(x0*y0*z1)*uj[3]+(x0*y1*z0)*uj[4]+(x1*y1*z0)*uj[5]+(x1*y1*z1)*uj[6]+(x0*y1*z1)*uj[7]; // perform trilinear interpolation
			c2 = shading(rainbow_color(255.0f*def_scale_u*length(ui)), p2+offset, normal, camera_cache); // rainbow_color(255.0f*def_scale_u*length(ui));
		}
		draw_triangle_interpolated(p0+offset, p1+offset, p2+offset, c0, c1, c2, camera_cache, bitmap, zbuffer); // draw triangle with interpolated colors
	}
//}

#endif // Graphics