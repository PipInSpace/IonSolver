
#if defined(FP16S) || defined(FP16C)
#define fpxx ushort
#else // FP32
#define fpxx float
#endif // FP32

//THESE ARE TO BE REMOVED
#define def_Nx 20u
#define def_Ny 20u
#define def_Nz 20u
#define def_N 8000ul

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
#define def_ke 8.9875517923E9f
#define def_charge 0.1f // Electric charge of a cell

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
#define ELECTRIC_FORCE

#define GRAPHICS
#define def_streamline_sparse 4u
#define def_streamline_length 128u
#define def_screen_width 1920u
#define def_screen_height 1080u
#define def_scale_u 1.0f
#define def_scale_Q_min 0.0001f
#define def_background_color 0x000000

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

// Helper functions 
float sq(const float x) {
	return x*x;
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
float cbmagnitude(uint3 v){
	return sq(v.x) + sq(v.y) + sq(v.z);
}

// Line3D OpenCL C version (c) Moritz Lehmann
//Graphics Helper functions:
// draw_point(...)    : draw 3D pixel
// draw_circle(...)   : draw 3D circle
// draw_line(...)     : draw 3D line
// draw_triangle(...) : draw 3D triangle
// iron_color(...)    : convert float in [0,255] to iron spectrum int color
// graphics_clear()   : kernel to reset bitmap and zbuffer
#ifdef GRAPHICS
int color_average(const int c1, const int c2) { // (c1+c2)/s
	const uchar4 cc1=as_uchar4(c1), cc2=as_uchar4(c2);
	return as_int((uchar4)((uchar)((cc1.x+cc2.x)/2u), (uchar)((cc1.y+cc2.y)/2u), (uchar)((cc1.z+cc2.z)/2u), (uchar)0u));
}
int color_mix_3(const int c0, const int c1, const int c2, const float w0, const float w1, const float w2) { // w1*c1+w2*c2+w3*c3, w0+w1+w2 = 1
	const uchar4 cc0=as_uchar4(c0), cc1=as_uchar4(c1), cc2=as_uchar4(c2);
	const float3 fc0=(float3)((float)cc0.x, (float)cc0.y, (float)cc0.z),  fc1=(float3)((float)cc1.x, (float)cc1.y, (float)cc1.z), fc2=(float3)((float)cc2.x, (float)cc2.y, (float)cc2.z);
	const float3 fcm = fma(w0, fc0, fma(w1, fc1, fma(w2, fc2, (float3)(0.5f, 0.5f, 0.5f))));
	return as_int((uchar4)((uchar)fcm.x, (uchar)fcm.y, (uchar)fcm.z, (uchar)0u));
}
int shading(const int c, const float3 p, const float3 normal, const float* camera_cache) {
    const float zoom = camera_cache[ 0]; // fetch camera parameters (rotation matrix, camera position, etc.)
	const float  dis = camera_cache[ 1];
	const float3 pos = (float3)(camera_cache[ 2], camera_cache[ 3], camera_cache[ 4])-(float3)(def_domain_offset_x, def_domain_offset_y, def_domain_offset_z);
	const float3 Rz  = (float3)(camera_cache[11], camera_cache[12], camera_cache[13]);
	const float3 d = p-Rz*(dis/zoom)-pos; // distance vector between p and camera position
    const float nl2 = sq(normal.x)+sq(normal.y)+sq(normal.z); // only one native_sqrt instead of two
    const float dl2 = sq(d.x)+sq(d.y)+sq(d.z);
	const float br = max(1.5f*fabs(dot(normal, d))*rsqrt(nl2*dl2), 0.3f);
	const float cr=(float)((c>>16)&255), cg=(float)((c>>8)&255), cb=(float)(c&255);
    return min((int)(br*cr), 255)<<16|min((int)(br*cg), 255)<<8|min((int)(br*cb), 255);
}
bool is_off_screen(const int x, const int y, const int stereo) {
	switch(stereo) {
		default: return x<                 0||x>=def_screen_width  ||y<0||y>=def_screen_height; // entire screen
		case -1: return x<                 0||x>=def_screen_width/2||y<0||y>=def_screen_height; // left half
		case +1: return x<def_screen_width/2||x>=def_screen_width  ||y<0||y>=def_screen_height; // right half
	}
}
void draw(const int x, const int y, const float z, const int color, global int* bitmap, volatile global int* zbuffer, const int stereo) {
	const int index=x+y*def_screen_width, iz=(int)(z*(2147483647.0f/10000.0f)); // use int z-buffer and atomic_max to minimize noise in image
	#ifndef GRAPHICS_TRANSPARENCY
		if(!is_off_screen(x, y, stereo)&&iz>atomic_max(&zbuffer[index], iz)) bitmap[index] = color; // only draw if point is on screen and first in zbuffer
	#else
		if(!is_off_screen(x, y, stereo)) { // transparent rendering (not quite order-independent transparency, but elegant solution for order-reversible transparency which is good enough here)
			const float transparency = GRAPHICS_TRANSPARENCY;
			const uchar4 cc4=as_uchar4(color), cb4=as_uchar4(def_background_color);
			const float3 fc = (float3)((float)cc4.x, (float)cc4.y, (float)cc4.z); // new pixel color that is behind topmost drawn pixel color
			const float3 fb = (float3)((float)cb4.x, (float)cb4.y, (float)cb4.z); // background color
			const bool is_front = iz>atomic_max(&zbuffer[index], iz);
			const uchar4 cp4 = as_uchar4(bitmap[index]);
			const float3 fp = (float3)((float)cp4.x, (float)cp4.y, (float)cp4.z); // current pixel color
			const int draw_count = (int)cp4.w; // use alpha color value to store how often the pixel has been over-drawn already
			const float3 fn = fp+(1.0f-transparency)*( is_front ? fc-fp : pown(transparency, draw_count)*(fc-fb)); // black magic: either over-draw colors back-to-front, or add back colors as correction terms
			bitmap[index] = as_int((uchar4)((uchar)clamp(fn.x+0.5f, 0.0f, 255.0f), (uchar)clamp(fn.y+0.5f, 0.0f, 255.0f), (uchar)clamp(fn.z+0.5f, 0.0f, 255.0f), (uchar)min(draw_count+1, 255)));
		}
	#endif
}
bool convert(int* rx, int* ry, float* rz, const float3 p, const float* camera_cache, const int stereo) { // 3D -> 2D
	const float zoom = camera_cache[0]; // fetch camera parameters (rotation matrix, camera position, etc.)
	const float  dis = camera_cache[1];
	const float3 pos = (float3)(camera_cache[ 2], camera_cache[ 3], camera_cache[ 4])-(float3)(def_domain_offset_x, def_domain_offset_y, def_domain_offset_z);
	const float3 Rx  = (float3)(camera_cache[ 5], camera_cache[ 6], camera_cache[ 7]);
	const float3 Ry  = (float3)(camera_cache[ 8], camera_cache[ 9], camera_cache[10]);
	const float3 Rz  = (float3)(camera_cache[11], camera_cache[12], camera_cache[13]);
	const float eye_distance = vload_half(28, (half*)camera_cache);
	float3 t, r;
	t = p-pos-((float)stereo*eye_distance/zoom)*(float3)(Rx.x, Rx.y, 0.0f); // transformation
	r.z = dot(Rz, t); // z-position for z-buffer
	const float rs = zoom*dis/(dis-r.z*zoom); // perspective (reciprocal is more efficient)
	if(rs<=0.0f) return false; // point is behins camera
	const float tv = ((as_int(camera_cache[14])>>30)&0x1)&&stereo!=0 ? 0.5f : 1.0f;
	r.x = (dot(Rx, t)*rs+(float)stereo*eye_distance)*tv+(0.5f+(float)stereo*0.25f)*(float)def_screen_width; // x position on screen
	r.y =  dot(Ry, t)*rs+0.5f*(float)def_screen_height; // y position on screen
	*rx = (int)(r.x+0.5f);
	*ry = (int)(r.y+0.5f);
	*rz = r.z;
	return true;
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
void convert_triangle(float3 p0, float3 p1, float3 p2, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer, const int stereo) { // 3D -> 2D
	int r0x, r0y, r1x, r1y, r2x, r2y; float r0z, r1z, r2z;
	if(convert(&r0x, &r0y, &r0z, p0, camera_cache, stereo) && convert(&r1x, &r1y, &r1z, p1, camera_cache, stereo) && convert(&r2x, &r2y, &r2z, p2, camera_cache, stereo)
		&& !(is_off_screen(r0x, r0y, stereo) && is_off_screen(r1x, r1y, stereo) && is_off_screen(r2x, r2y, stereo))) { // cancel drawing if all points are off screen
		if(r0x*(r1y-r2y)+r1x*(r2y-r0y)+r2x*(r0y-r1y)>40000 || (r0y==r1y&&r0y==r2y)) return; // return for large triangle area or degenerate triangles
		//if(r1x*r0y+r2x*r1y+r0x*r2y>=r0x*r1y+r1x*r2y+r2x*r0y) return; // clockwise backface culling
		if(r0y>r1y) { const int xt = r0x; const int yt = r0y; r0x = r1x; r0y = r1y; r1x = xt; r1y = yt; } // sort vertices ascending by y
		if(r0y>r2y) { const int xt = r0x; const int yt = r0y; r0x = r2x; r0y = r2y; r2x = xt; r2y = yt; }
		if(r1y>r2y) { const int xt = r1x; const int yt = r1y; r1x = r2x; r1y = r2y; r2x = xt; r2y = yt; }
		const float z = (r0z+r1z+r2z)/3.0f; // approximate triangle z position for each pixel to be equal
		for(int y=r0y; y<r1y; y++) { // Bresenham algorithm (lower triangle half)
			const int xA = r0x+(r2x-r0x)*(y-r0y)/(r2y-r0y);
			const int xB = r0x+(r1x-r0x)*(y-r0y)/(r1y-r0y);
			for(int x=min(xA, xB); x<max(xA, xB); x++) {
				draw(x, y, z, color, bitmap, zbuffer, stereo);
			}
		}
		for(int y=r1y; y<r2y; y++) { // Bresenham algorithm (upper triangle half)
			const int xA = r0x+(r2x-r0x)*(y-r0y)/(r2y-r0y);
			const int xB = r1x+(r2x-r1x)*(y-r1y)/(r2y-r1y);
			for(int x=min(xA, xB); x<max(xA, xB); x++) {
				draw(x, y, z, color, bitmap, zbuffer, stereo);
			}
		}
	}
}
void convert_triangle_interpolated(float3 p0, float3 p1, float3 p2, int c0, int c1, int c2, const float* camera_cache, global int* bitmap, global int* zbuffer, const int stereo) { // 3D -> 2D
	int r0x, r0y, r1x, r1y, r2x, r2y; float r0z, r1z, r2z;
	if(convert(&r0x, &r0y, &r0z, p0, camera_cache, stereo) && convert(&r1x, &r1y, &r1z, p1, camera_cache, stereo) && convert(&r2x, &r2y, &r2z, p2, camera_cache, stereo)
		&& !(is_off_screen(r0x, r0y, stereo) && is_off_screen(r1x, r1y, stereo) && is_off_screen(r2x, r2y, stereo))) { // cancel drawing if all points are off screen
		if(r0x*(r1y-r2y)+r1x*(r2y-r0y)+r2x*(r0y-r1y)>40000 || (r0y==r1y&&r0y==r2y)) return; // return for large triangle area or degenerate triangles
		//if(r1x*r0y+r2x*r1y+r0x*r2y>=r0x*r1y+r1x*r2y+r2x*r0y) return; // clockwise backface culling
		if(r0y>r1y) { const int xt = r0x; const int yt = r0y; r0x = r1x; r0y = r1y; r1x = xt; r1y = yt; const int ct = c0; c0 = c1; c1 = ct; } // sort vertices ascending by y
		if(r0y>r2y) { const int xt = r0x; const int yt = r0y; r0x = r2x; r0y = r2y; r2x = xt; r2y = yt; const int ct = c0; c0 = c2; c2 = ct; }
		if(r1y>r2y) { const int xt = r1x; const int yt = r1y; r1x = r2x; r1y = r2y; r2x = xt; r2y = yt; const int ct = c1; c1 = c2; c2 = ct; }
		const float z = (r0z+r1z+r2z)/3.0f; // approximate triangle z position for each pixel to be equal
		const float d = (float)((r1y-r2y)*(r0x-r2x)+(r2x-r1x)*(r0y-r2y));
		for(int y=r0y; y<r1y; y++) { // Bresenham algorithm (lower triangle half)
			const int xA = r0x+(r2x-r0x)*(y-r0y)/(r2y-r0y);
			const int xB = r0x+(r1x-r0x)*(y-r0y)/(r1y-r0y);
			for(int x=min(xA, xB); x<max(xA, xB); x++) {
				const float w0 = (float)((r1y-r2y)*(x-r2x)+(r2x-r1x)*(y-r2y))/d; // barycentric coordinates
				const float w1 = (float)((r2y-r0y)*(x-r2x)+(r0x-r2x)*(y-r2y))/d;
				const float w2 = 1.0f-w0-w1;
				const int color = color_mix_3(c0, c1, c2, w0, w1, w2); // interpolate color
				draw(x, y, z, color, bitmap, zbuffer, stereo);
			}
		}
		for(int y=r1y; y<r2y; y++) { // Bresenham algorithm (upper triangle half)
			const int xA = r0x+(r2x-r0x)*(y-r0y)/(r2y-r0y);
			const int xB = r1x+(r2x-r1x)*(y-r1y)/(r2y-r1y);
			for(int x=min(xA, xB); x<max(xA, xB); x++) {
				const float w0 = (float)((r1y-r2y)*(x-r2x)+(r2x-r1x)*(y-r2y))/d; // barycentric coordinates
				const float w1 = (float)((r2y-r0y)*(x-r2x)+(r0x-r2x)*(y-r2y))/d;
				const float w2 = 1.0f-w0-w1;
				const int color = color_mix_3(c0, c1, c2, w0, w1, w2); // interpolate color
				draw(x, y, z, color, bitmap, zbuffer, stereo);
			}
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
void draw_triangle(const float3 p0, const float3 p1, const float3 p2, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer) { // 3D -> 2D
	const bool vr = (as_int(camera_cache[14])>>31)&0x1;
	if(!vr) {
		convert_triangle(p0, p1, p2, color, camera_cache, bitmap, zbuffer,  0);
	} else {
		convert_triangle(p0, p1, p2, color, camera_cache, bitmap, zbuffer, -1); // left eye
		convert_triangle(p0, p1, p2, color, camera_cache, bitmap, zbuffer, +1); // right eye
	}
}
void draw_triangle_interpolated(const float3 p0, const float3 p1, const float3 p2, const int c0, const int c1, const int c2, const float* camera_cache, global int* bitmap, global int* zbuffer) { // 3D -> 2D
	const bool vr = (as_int(camera_cache[14])>>31)&0x1;
	if(!vr) {
		convert_triangle_interpolated(p0, p1, p2, c0, c1, c2, camera_cache, bitmap, zbuffer,  0);
	} else {
		convert_triangle_interpolated(p0, p1, p2, c0, c1, c2, camera_cache, bitmap, zbuffer, -1); // left eye
		convert_triangle_interpolated(p0, p1, p2, c0, c1, c2, camera_cache, bitmap, zbuffer, +1); // right eye
	}
}
kernel void graphics_clear(global int* bitmap, global int* zbuffer) {
	const uint n = get_global_id(0);
	bitmap[n] = def_background_color; // black background = 0x000000, use 0xFFFFFF for white background
	zbuffer[n] = -2147483648;
}
constant ushort edge_table_data[128] = { // source: Paul Bourke, http://paulbourke.net/geometry/polygonise/, mirror symmetry applied, makes marching-cubes 31% faster
	0x000, 0x109, 0x203, 0x30A, 0x406, 0x50F, 0x605, 0x70C, 0x80C, 0x905, 0xA0F, 0xB06, 0xC0A, 0xD03, 0xE09, 0xF00,
	0x190, 0x099, 0x393, 0x29A, 0x596, 0x49F, 0x795, 0x69C, 0x99C, 0x895, 0xB9F, 0xA96, 0xD9A, 0xC93, 0xF99, 0xE90,
	0x230, 0x339, 0x033, 0x13A, 0x636, 0x73F, 0x435, 0x53C, 0xA3C, 0xB35, 0x83F, 0x936, 0xE3A, 0xF33, 0xC39, 0xD30,
	0x3A0, 0x2A9, 0x1A3, 0x0AA, 0x7A6, 0x6AF, 0x5A5, 0x4AC, 0xBAC, 0xAA5, 0x9AF, 0x8A6, 0xFAA, 0xEA3, 0xDA9, 0xCA0,
	0x460, 0x569, 0x663, 0x76A, 0x066, 0x16F, 0x265, 0x36C, 0xC6C, 0xD65, 0xE6F, 0xF66, 0x86A, 0x963, 0xA69, 0xB60,
	0x5F0, 0x4F9, 0x7F3, 0x6FA, 0x1F6, 0x0FF, 0x3F5, 0x2FC, 0xDFC, 0xCF5, 0xFFF, 0xEF6, 0x9FA, 0x8F3, 0xBF9, 0xAF0,
	0x650, 0x759, 0x453, 0x55A, 0x256, 0x35F, 0x055, 0x15C, 0xE5C, 0xF55, 0xC5F, 0xD56, 0xA5A, 0xB53, 0x859, 0x950,
	0x7C0, 0x6C9, 0x5C3, 0x4CA, 0x3C6, 0x2CF, 0x1C5, 0x0CC, 0xFCC, 0xEC5, 0xDCF, 0xCC6, 0xBCA, 0xAC3, 0x9C9, 0x8C0
};
constant uchar triangle_table_data[1920] = { // source: Paul Bourke, http://paulbourke.net/geometry/polygonise/, termination value 15, bit packed
	255,255,255,255,255,255,255, 15, 56,255,255,255,255,255,255, 16,249,255,255,255,255,255, 31, 56,137,241,255,255,255,255, 33,250,255,255,255,255,255, 15, 56, 33,250,255,255,255,255, 41, 10,146,
	255,255,255,255, 47, 56,162,168,137,255,255,255,179,242,255,255,255,255,255, 15, 43,184,240,255,255,255,255,145, 32,179,255,255,255,255, 31, 43,145,155,184,255,255,255,163,177, 58,255,255,255,
	255, 15, 26,128,138,171,255,255,255,147, 48,155,171,249,255,255,159,168,138,251,255,255,255,255,116,248,255,255,255,255,255, 79,  3, 55,244,255,255,255,255, 16,137,116,255,255,255,255, 79,145,
	116,113, 19,255,255,255, 33,138,116,255,255,255,255, 63,116,  3, 20,162,255,255,255, 41,154, 32, 72,247,255,255, 47,154,146, 39, 55,151,244,255, 72, 55, 43,255,255,255,255,191,116, 43, 36, 64,
	255,255,255,  9,129,116, 50,251,255,255, 79,183, 73,155, 43, 41,241,255,163, 49,171,135,244,255,255, 31,171, 65, 27, 64,183,244,255,116,152,176,185,186, 48,255, 79,183,180,153,171,255,255,255,
	 89,244,255,255,255,255,255,159, 69,128,243,255,255,255,255, 80, 20,  5,255,255,255,255,143, 69, 56, 53, 81,255,255,255, 33,154, 69,255,255,255,255, 63,128, 33, 74, 89,255,255,255, 37, 90, 36,
	  4,242,255,255, 47, 90, 35, 53, 69, 67,248,255, 89, 36,179,255,255,255,255, 15, 43,128, 75, 89,255,255,255, 80,  4, 81, 50,251,255,255, 47, 81, 82, 40,184,132,245,255, 58,171, 49, 89,244,255,
	255, 79, 89,128,129, 26,184,250,255, 69, 80,176,181,186, 48,255, 95,132,133,170,184,255,255,255,121, 88,151,255,255,255,255,159,  3, 89, 83, 55,255,255,255,112,  8,113, 81,247,255,255, 31, 53,
	 83,247,255,255,255,255,121,152,117, 26,242,255,255,175, 33, 89, 80,  3,117,243,255,  8,130, 82, 88,167, 37,255, 47, 90, 82, 51,117,255,255,255,151,117,152,179,242,255,255,159,117,121,146,  2,
	114,251,255, 50, 11,129,113, 24,117,255,191, 18, 27,119, 81,255,255,255, 89,136,117, 26,163,179,255, 95,  7,  5,121, 11,  1,186, 10,171,176, 48, 90,128,112,117,176, 90,183,245,255,255,255,255,
	106,245,255,255,255,255,255, 15, 56,165,246,255,255,255,255,  9, 81,106,255,255,255,255, 31, 56,145, 88,106,255,255,255, 97, 37, 22,255,255,255,255, 31, 86, 33, 54,128,255,255,255,105,149, 96,
	 32,246,255,255, 95,137,133, 82, 98, 35,248,255, 50,171, 86,255,255,255,255,191,128, 43,160, 86,255,255,255, 16, 41,179,165,246,255,255, 95,106,145,146, 43,137,251,255, 54,107, 53, 21,243,255,
	255, 15,184,176,  5, 21,181,246,255,179,  6, 99, 96,  5,149,255,111,149,150,187,137,255,255,255,165, 70,135,255,255,255,255, 79,  3,116, 99,165,255,255,255,145, 80,106, 72,247,255,255,175, 86,
	145, 23, 55,151,244,255, 22, 98, 21,116,248,255,255, 31, 82, 37, 54, 64, 67,247,255, 72,151, 80, 96,  5, 98,255,127,147,151, 52,146,149, 38,150,179,114, 72,106,245,255,255, 95,106,116, 66,  2,
	114,251,255, 16, 73,135, 50, 91,106,255,159, 18,185,146,180,183, 84,106, 72, 55, 91, 83, 81,107,255, 95,177,181, 22,176,183,  4,180, 80,  9, 86, 48,182, 54, 72,103,149,150, 75,151,183,249,255,
	 74,105,164,255,255,255,255, 79,106,148, 10, 56,255,255,255, 10,161,  6, 70,240,255,255,143, 19, 24,134, 70, 22,250,255, 65, 25, 66, 98,244,255,255, 63,128, 33, 41,148, 98,244,255, 32, 68, 98,
	255,255,255,255,143, 35, 40, 68, 98,255,255,255, 74,169, 70, 43,243,255,255, 15, 40,130, 75,169,164,246,255,179,  2, 97, 96,100,161,255,111, 20, 22, 74, 24, 18,139, 27,105,148, 99, 25,179, 54,
	255,143, 27, 24,176, 22, 25,100, 20,179, 54,  6, 96,244,255,255,111,132,107,248,255,255,255,255,167,118,168,152,250,255,255, 15, 55,160,  7,169,118,250,255,106, 23,122,113, 24,  8,255,175,118,
	122, 17, 55,255,255,255, 33, 22,134,129,137,118,255, 47,150,146, 97,151,144,115,147,135,112, 96,  6,242,255,255,127, 35,118,242,255,255,255,255, 50,171,134,138,137,118,255, 47,112,114, 11,121,
	118,154,122,129, 16,135,161,103,167, 50,187, 18, 27,167, 22,118,241,255,152,134,118, 25,182, 54, 49,  6, 25,107,247,255,255,255,255,135,112, 96,179,176,  6,255,127,107,255,255,255,255,255,255,
	103,251,255,255,255,255,255, 63,128,123,246,255,255,255,255, 16,185,103,255,255,255,255,143,145, 56,177,103,255,255,255, 26, 98,123,255,255,255,255, 31,162,  3,104,123,255,255,255,146, 32,154,
	182,247,255,255,111,123,162,163, 56,154,248,255, 39, 99,114,255,255,255,255,127,128,103, 96,  2,255,255,255,114, 38,115, 16,249,255,255, 31, 38,129, 22,137,120,246,255,122,166,113, 49,247,255,
	255,175,103,113, 26,120,  1,248,255, 48,  7,167,160,105,122,255,127,166,167,136,154,255,255,255,134,180,104,255,255,255,255, 63,182,  3,  6,100,255,255,255,104,139,100,  9,241,255,255,159,100,
	105,147, 19, 59,246,255,134,100,139,162,241,255,255, 31,162,  3, 11,182, 64,246,255,180, 72,182, 32, 41,154,255,175, 57, 58,146, 52, 59, 70, 54, 40,131, 36,100,242,255,255, 15, 36,100,242,255,
	255,255,255,145, 32, 67, 66, 70,131,255, 31, 73, 65, 34,100,255,255,255, 24,131, 22, 72,102, 26,255,175,  1, 10,102, 64,255,255,255,100, 67,131,166,  3,147,154,163, 73,166,244,255,255,255,255,
	148,117,182,255,255,255,255, 15, 56,148,181,103,255,255,255,  5, 81,  4,103,251,255,255,191,103, 56, 52, 69, 19,245,255, 89,164, 33,103,251,255,255,111,123, 33, 10, 56,148,245,255,103, 91,164,
	 36, 74, 32,255, 63,132, 83, 52, 82, 90,178,103, 39,115, 38, 69,249,255,255,159, 69,128,  6, 38,134,247,255, 99, 50,103, 81, 80,  4,255,111,130,134, 39,129,132, 21,133, 89,164, 97,113, 22,115,
	255, 31,166,113, 22,112,120,144, 69,  4, 74, 90, 48,106,122,115,122,166,167, 88,164,132,250,255,150,101,155,139,249,255,255, 63,182, 96,  3,101,144,245,255,176,  8,181, 16, 85,182,255,111, 59,
	 54, 85, 19,255,255,255, 33,154,181,185,184,101,255, 15, 59, 96, 11,105,101, 25,162,139,181,101,  8,165, 37, 32,101, 59, 54, 37, 58, 90,243,255,133, 89,130,101, 50, 40,255,159,101,105,  0, 38,
	255,255,255, 81, 24,  8,101, 56, 40, 38, 24,101, 18,246,255,255,255,255, 49, 22,166,131, 86,150,152,166,  1, 10,150,  5,101,240,255, 48, 88,166,255,255,255,255,175,101,255,255,255,255,255,255,
	 91,122,181,255,255,255,255,191,165,123,133,  3,255,255,255,181, 87,186,145,240,255,255,175, 87,186,151, 24, 56,241,255, 27,178, 23, 87,241,255,255, 15, 56, 33, 23, 87, 39,251,255,121,149,114,
	  9, 34,123,255,127, 37, 39, 91, 41, 35,152, 40, 82, 42, 83,115,245,255,255,143,  2, 88,130, 87, 42,245,255,  9, 81, 58, 53, 55, 42,255,159, 40, 41,129, 39, 42,117, 37, 49, 53, 87,255,255,255,
	255, 15,120,112, 17, 87,255,255,255,  9,147, 83, 53,247,255,255,159,120,149,247,255,255,255,255,133, 84,138,186,248,255,255, 95, 64,181, 80,186, 59,240,255, 16,137,164,168,171, 84,255,175, 75,
	 74,181, 67, 73, 49, 65, 82, 33, 88,178, 72,133,255, 15,180,176, 67,181,178, 81,177, 32,  5,149,178, 69,133,139,149, 84,178,243,255,255,255,255, 82, 58, 37, 67, 53, 72,255, 95, 42, 37, 68,  2,
	255,255,255,163, 50,165,131, 69,133, 16, 89, 42, 37, 20, 41, 73,242,255, 72,133, 53, 83,241,255,255, 15, 84,  1,245,255,255,255,255, 72,133, 53,  9,  5, 83,255,159, 84,255,255,255,255,255,255,
	180, 71,185,169,251,255,255, 15, 56,148,151,123,169,251,255,161, 27, 75, 65,112,180,255, 63, 65, 67, 24, 74, 71,171, 75,180,151, 75, 41,155, 33,255,159, 71,185,151,177,178,  1, 56,123,180, 36,
	 66,240,255,255,191, 71, 75,130, 67, 35,244,255,146, 42,151, 50,119,148,255,159,122,121,164,114,120, 32,112,115, 58, 42, 71, 26, 10,  4, 26, 42,120,244,255,255,255,255,148, 65,113, 23,243,255,
	255, 79, 25, 20,  7, 24,120,241,255,  4,115, 52,255,255,255,255, 79,120,255,255,255,255,255,255,169,168,139,255,255,255,255, 63,144,147,187,169,255,255,255, 16, 10,138,168,251,255,255, 63,161,
	 59,250,255,255,255,255, 33, 27,155,185,248,255,255, 63,144,147, 27,146,178,249,255, 32,139,176,255,255,255,255, 63,178,255,255,255,255,255,255, 50, 40,168,138,249,255,255,159, 42,144,242,255,
	255,255,255, 50, 40,168, 16, 24,138,255, 31, 42,255,255,255,255,255,255, 49,152,129,255,255,255,255, 15, 25,255,255,255,255,255,255, 48,248,255,255,255,255,255,255,255,255,255,255,255,255,255
};
ushort edge_table(const uint i) {
	return edge_table_data[i<128u?i:255u-i];
}
uchar triangle_table(const uint i) {
	return (triangle_table_data[i/2u]>>(4u*(i%2u)))&0xF;
}
float3 interpolate_vertex(const float3 p1, const float3 p2, const float v1, const float v2, const float iso) { // linearly interpolate position where isosurface cuts an edge between 2 vertices
	const float w = (iso-v1)/(v2-v1);
	return (1.0f-w)*p1+w*p2;
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
bool is_above_plane(const float3 point, const float3 plane_p, const float3 plane_n) {
	return dot(point-plane_p, plane_n)>=0.0f;
}
bool is_in_camera_frustrum(const float3 p, const float* camera_cache) { // returns true if point is located in camera frustrum
	const float zoom = camera_cache[0]; // fetch camera parameters (rotation matrix, camera position, etc.)
	const float  dis = camera_cache[1];
	const float3 pos = (float3)(camera_cache[ 2], camera_cache[ 3], camera_cache[ 4])-(float3)(def_domain_offset_x, def_domain_offset_y, def_domain_offset_z);
	const float3 Rx  = (float3)(camera_cache[ 5], camera_cache[ 6], camera_cache[ 7]);
	const float3 Ry  = (float3)(camera_cache[ 8], camera_cache[ 9], camera_cache[10]);
	const float3 Rz  = (float3)(camera_cache[11], camera_cache[12], camera_cache[13]);
	const bool   vr  = (as_int(camera_cache[14])>>31)&0x1;
	const float  rtv = (as_int(camera_cache[14])>>30)&0x1 ? 2.0f : 1.0f;
	const float3 p0 = (float3)(0.0f, 0.0f, dis/zoom);
	const float3 camera_center = Rx*p0.x+Ry*p0.y+Rz*p0.z+pos; // reverse rotation and reverse transformation of p0
	const float x_left   = !vr ? (float)(-(int)def_screen_width/2  ) : ((float)(-(int)def_screen_width/2  )+(float)(def_screen_width/4u))*rtv;
	const float x_right  = !vr ? (float)( (int)def_screen_width/2-1) : ((float)( (int)def_screen_width/2-1)-(float)(def_screen_width/4u))*rtv;
	const float y_top    = (float)(-(int)def_screen_height/2 );
	const float y_bottom = (float)((int)def_screen_height/2-1);
	float3 r00 = p0+normalize((float3)(x_left , y_top   , -dis)); // get 4 edge vectors of frustrum, get_camray(...) inlined and redundant parts eliminated
	float3 r01 = p0+normalize((float3)(x_right, y_top   , -dis));
	float3 r10 = p0+normalize((float3)(x_left , y_bottom, -dis));
	float3 r11 = p0+normalize((float3)(x_right, y_bottom, -dis));
	r00 = Rx*r00.x+Ry*r00.y+Rz*r00.z+pos-camera_center; // reverse rotation and reverse transformation of r00
	r01 = Rx*r01.x+Ry*r01.y+Rz*r01.z+pos-camera_center; // reverse rotation and reverse transformation of r01
	r10 = Rx*r10.x+Ry*r10.y+Rz*r10.z+pos-camera_center; // reverse rotation and reverse transformation of r10
	r11 = Rx*r11.x+Ry*r11.y+Rz*r11.z+pos-camera_center; // reverse rotation and reverse transformation of r11
	const float3 plane_n_top    = cross(r00, r01); // get 4 frustrum planes
	const float3 plane_n_bottom = cross(r11, r10);
	const float3 plane_n_left   = cross(r10, r00);
	const float3 plane_n_right  = cross(r01, r11);
	const float3 plane_p_top    = camera_center-2.0f*plane_n_top; // move frustrum planes outward by 2 units
	const float3 plane_p_bottom = camera_center-2.0f*plane_n_bottom;
	const float3 plane_p_left   = camera_center-(2.0f+8.0f*(float)vr)*plane_n_left; // move frustrum planes outward by 2 units, for stereoscopic rendering a bit more
	const float3 plane_p_right  = camera_center-(2.0f+8.0f*(float)vr)*plane_n_right;
	return is_above_plane(p, plane_p_top, plane_n_top)&&is_above_plane(p, plane_p_bottom, plane_n_bottom)&&is_above_plane(p, plane_p_left, plane_n_left)&&is_above_plane(p, plane_p_right, plane_n_right);
}

#endif

float3 position(const uint3 xyz) { // 3D coordinates to 3D position
	return (float3)((float)xyz.x+0.5f-0.5f*(float)def_Nx, (float)xyz.y+0.5f-0.5f*(float)def_Ny, (float)xyz.z+0.5f-0.5f*(float)def_Nz);
}

uint3 coordinates(const uint n) { // disassemble 1D index to 3D coordinates (n -> x,y,z)
	const uint t = n%(def_Nx*def_Ny);
	return (uint3)(t%def_Nx, t/def_Nx, n/(def_Nx*def_Ny)); // n = x+(y+z*Ny)*Nx
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

#ifdef ELECTRIC_FORCE
// we need to optimize this
// n: cell id
// q: float array for charges
// E: electric field at n
void calculate_E(const uint n, const global float* q, global float* E) {// uses coulomb's law https://en.wikipedia.org/wiki/Coulomb%27s_law
	const float3 c_n = convert_float3(coordinates(n));
	float3 e_i;
	for(uint i = 0; i < def_N; i++){
		const float3 c_i = convert_float3(coordinates(i));
		//const float3 e_i = def_ke * convert_float3((coord_n - coord_i)) / cbmagnitude(coord_n - coord_i); // coulomb's law
		e_i += q[i] / (sq(c_n.x-c_i.x)+sq(c_n.y-c_i.y)+sq(c_n.z-c_i.z)) * fast_normalize(c_n - c_i);
	}
	E[n					] = e_i.x * def_ke;
	E[(ulong)n+def_N	] = e_i.y * def_ke;
	E[(ulong)n+def_N*2ul] = e_i.z * def_ke;
}
#endif

__kernel void stream_collide(global fpxx* fi, global float* rho, global float* u, global uchar* flags, const ulong t, const float fx, const float fy, const float fz 
#ifdef FORCE_FIELD
, const global float* F 
#endif // FORCE_FIELD
#ifdef ELECTRIC_FORCE
, global float* E
#endif // ELECTRIC_FORCE
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

	#ifdef ELECTRIC_FORCE
	//{ // separate block to avoid variable name conflicts
	//	fxn += E[                 n] * rhon * def_charge; // apply electric field * charge = force
	//	fyn += E[    def_N+(ulong)n] * rhon * def_charge;
	//	fzn += E[2ul*def_N+(ulong)n] * rhon * def_charge;
	//}
	#endif// ELECTRIC_FORCE

    float Fin[def_velocity_set]; // forcing terms

	#ifdef FORCE_FIELD
	{ // separate block to avoid variable name conflicts
		fxn += F[                 n]; // apply force field
		fyn += F[    def_N+(ulong)n];
		fzn += F[2ul*def_N+(ulong)n];
	}
	#endif

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
		rho[               n] = rhon; // update density field
		u[                 n] = uxn; // update velocity field
		u[    def_N+(ulong)n] = uyn;
		u[2ul*def_N+(ulong)n] = uzn;
	#endif // UPDATE_FIELDS
	#else // EQUILIBRIUM_BOUNDARIES
	#ifdef UPDATE_FIELDS
		if(flagsn_bo!=TYPE_E) { // only update fields for non-TYPE_E cells
			rho[               n] = rhon; // update density field
			u[                 n] = uxn; // update velocity field
			u[    def_N+(ulong)n] = uyn;
			u[2ul*def_N+(ulong)n] = uzn;
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
#ifdef ELECTRIC_FORCE
, global float* q, global float* E
#endif
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
	#ifdef ELECTRIC_FORCE
		//calculate_E(n, q, E);
	#endif // ELECTRIC FORCE
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
	// update charge field (later)


}

// Graphics code
#ifdef GRAPHICS
int iron_color(float x) { // coloring scheme (float 0-255 -> int color)
	x = clamp(360.0f-x*360.0f/255.0f, 0.0f, 360.0f);
	float r=255.0f, g=0.0f, b=0.0f;
	if(x<60.0f) { // white - yellow
		g = 255.0f;
		b = 255.0f-255.0f*x/60.0f;
	} else if(x<180.0f) { // yellow - red
		g = 255.0f-255.0f*(x-60.0f)/120.0f;
	} else if(x<270.0f) { // red - violet
		r = 255.0f-255.0f*(x-180.0f)/180.0f;
		b = 255.0f*(x-180.0f)/90.0f;
	} else { // violet - black
		r = 255.0f-255.0f*(x-180.0f)/180.0f;
		b = 255.0f-255.0f*(x-270.0f)/90.0f;
	}
	return (((int)r)<<16)|(((int)g)<<8)|((int)b);
}
int rainbow_color(float x) { // coloring scheme (float 0-255 -> int color)
	x = clamp(360.0f-x*360.0f/255.0f, 0.0f, 360.0f);
	float r=0.0f, g=0.0f, b=0.0f; // black
	if(x<60.0f) { // red - yellow
		r = 255.0f;
		g = 255.0f*x/60.0f;
	} else if(x>=60.0f&&x<120.0f) { // yellow - green
		r = 255.0f-255.0f*(x-60.0f)/60.0f;
		g = 255.0f;
	} else if(x>=120.0f&&x<180.0f) { // green - cyan
		g = 255.0f;
		b = 255.0f*(x-120.0f)/60.0f;
	} else if(x>=180.0f&&x<240.0f) { // cyan - blue
		g = 255.0f-255.0f*(x-180.0f)/60.0f;
		b = 255.0f;
	} else if(x>=240.0f&&x<300.0f) { // blue - violet
		r = (255.0f*(x-240.0f)/60.0f)/2.0f;
		b = 255.0f;
	} else { // violet - black
		r = (255.0f-255.0f*(x-300.0f)/60.0f)/2.0f;
		b = 255.0f-255.0f*(x-300.0f)/60.0f;
	}
	return (((int)r)<<16)|(((int)g)<<8)|((int)b);
}
kernel void graphics_flags(const global uchar* flags, const global float* camera, global int* bitmap, global int* zbuffer) {
    const uint n = get_global_id(0);
    if(n>=(uint)def_N||is_halo(n)) return; // don't execute graphics_flags() on halo
    const uchar flagsn = flags[n]; // cache flags
    const uchar flagsn_bo = flagsn&TYPE_BO; // extract boundary flags
    if(flagsn==0u||flagsn==TYPE_G) return; // don't draw regular fluid cells
    //if(flagsn&TYPE_SU) return; // don't draw surface
    float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
    for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
    const uint3 xyz = coordinates(n);
    const float3 p = position(xyz);
	if(!is_in_camera_frustrum(p, camera_cache)) return; // skip loading LBM data if grid cell is not visible
	uint x0, xp, xm, y0, yp, ym, z0, zp, zm;
	calculate_indices(n, &x0, &xp, &xm, &y0, &yp, &ym, &z0, &zp, &zm);
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
	float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
	for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
	const float3 p = position(xyz);
	if(!is_in_camera_frustrum(p, camera_cache)) return; // skip loading LBM data if grid cell is not visible
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
    for(uint i=0u; i<tn; i++) { // TODO: make compatible with FORCE_FIELD
	    const float3 p0 = triangles[3u*i   ];
	    const float3 p1 = triangles[3u*i+1u];
	    const float3 p2 = triangles[3u*i+2u];
	    const float3 normal = normalize(cross(p1-p0, p2-p0));

        const int c = shading(0xDFDFDF, p+(p0+p1+p2)/3.0f, normal, camera_cache); // 0xDFDFDF;
		draw_triangle(p+p0, p+p1, p+p2, c, camera_cache, bitmap, zbuffer);
    }
}

kernel void graphics_field(const global uchar* flags, const global float* u, const global float* camera, global int* bitmap, global int* zbuffer, const int slice_mode, const int slice_x, const int slice_y, const int slice_z) {
	const uint n = get_global_id(0);
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute graphics_field() on halo
	const uint3 xyz = coordinates(n);
	const bool rx=(int)xyz.x!=slice_x, ry=(int)xyz.y!=slice_y, rz=(int)xyz.z!=slice_z;
	if((slice_mode==1&&rx)||(slice_mode==2&&ry)||(slice_mode==3&&rz)||(slice_mode==4&&rx&&rz)||(slice_mode==5&&rx&&ry&&rz)||(slice_mode==6&&ry&&rz)||(slice_mode==7&&rx&&ry)) return;
	float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
	for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
	const float3 p = position(xyz);
	if(!is_in_camera_frustrum(p, camera_cache)) return; // skip loading LBM data if grid cell is not visible

    if(flags[n]&(TYPE_S|TYPE_E|TYPE_I|TYPE_G)) return;

    float3 un = load_u(n, u); // cache velocity
    const float ul = length(un);
    if(def_scale_u*ul<0.1f) return; // don't draw lattice points where the velocity is lower than this threshold
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
    #else // D2Q9
        if(n>=(def_Nx/def_streamline_sparse)*(def_Ny/def_streamline_sparse)) return;
        const uint y = n/(def_Nx/def_streamline_sparse); // disassemble 1D index to 3D coordinates
        const uint x = n%(def_Nx/def_streamline_sparse);
        float3 p = ((float3)((float)def_streamline_sparse*((float)x+0.5f), (float)def_streamline_sparse*((float)y+0.5f), 0.5f))-0.5f*((float3)((float)def_Nx, (float)def_Ny, (float)def_Nz));
        const bool rx=fabs(p.x-slice.x)>0.5f*(float)def_streamline_sparse, ry=fabs(p.y-slice.y)>0.5f*(float)def_streamline_sparse, rz=true;
    #endif // D2Q9
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
	const float3 p = position(coordinates(n));
	float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
	for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
	if(!is_in_camera_frustrum(p, camera_cache)) return; // skip loading LBM data if grid cell is not visible
    float3 un = load_u(n, u); // cache velocity
    const float ul = length(un);
    const float Q = calculate_Q(n, u);
    if(Q<def_scale_Q_min||ul==0.0f) return; // don't draw lattice points where the velocity is very low
    const int c = rainbow_color(255.0f*def_scale_u*ul); // coloring by velocity
    draw_line(p-(0.5f/ul)*un, p+(0.5f/ul)*un, c, camera_cache, bitmap, zbuffer);
}

kernel void graphics_q(const global uchar* flags, const global float* u, const global float* camera, global int* bitmap, global int* zbuffer) {
    const uint n = get_global_id(0);
	const uint3 xyz = coordinates(n);
	if(xyz.x>=def_Nx-1u||xyz.y>=def_Ny-1u||xyz.z>=def_Nz-1u||is_halo_q(xyz)) return; // don't execute graphics_q_field() on marching-cubes halo
	const float3 p = position(xyz);
	float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
	for(uint i=0u; i<15u; i++) camera_cache[i] = camera[i];
	if(!is_in_camera_frustrum(p, camera_cache)) return; // skip loading LBM data if grid cell is not visible
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
	for(uint i=0u; i<tn; i++) {
		const float3 p0 = triangles[3u*i   ]; // triangle coordinates in [0,1] (local cell)
		const float3 p1 = triangles[3u*i+1u];
		const float3 p2 = triangles[3u*i+2u];
		const float3 normal = normalize(cross(p1-p0, p2-p0));
		int c0, c1, c2; {
			const float x1=p0.x, y1=p0.y, z1=p0.z, x0=1.0f-x1, y0=1.0f-y1, z0=1.0f-z1; // calculate interpolation factors
			const float3 ui = (x0*y0*z0)*uj[0]+(x1*y0*z0)*uj[1]+(x1*y0*z1)*uj[2]+(x0*y0*z1)*uj[3]+(x0*y1*z0)*uj[4]+(x1*y1*z0)*uj[5]+(x1*y1*z1)*uj[6]+(x0*y1*z1)*uj[7]; // perform trilinear interpolation
			c0 = shading(rainbow_color(255.0f*def_scale_u*length(ui)), p+p0, normal, camera_cache); // rainbow_color(255.0f*def_scale_u*length(ui));
		} {
			const float x1=p1.x, y1=p1.y, z1=p1.z, x0=1.0f-x1, y0=1.0f-y1, z0=1.0f-z1; // calculate interpolation factors
			const float3 ui = (x0*y0*z0)*uj[0]+(x1*y0*z0)*uj[1]+(x1*y0*z1)*uj[2]+(x0*y0*z1)*uj[3]+(x0*y1*z0)*uj[4]+(x1*y1*z0)*uj[5]+(x1*y1*z1)*uj[6]+(x0*y1*z1)*uj[7]; // perform trilinear interpolation
			c1 = shading(rainbow_color(255.0f*def_scale_u*length(ui)), p+p1, normal, camera_cache); // rainbow_color(255.0f*def_scale_u*length(ui));
		} {
			const float x1=p2.x, y1=p2.y, z1=p2.z, x0=1.0f-x1, y0=1.0f-y1, z0=1.0f-z1; // calculate interpolation factors
			const float3 ui = (x0*y0*z0)*uj[0]+(x1*y0*z0)*uj[1]+(x1*y0*z1)*uj[2]+(x0*y0*z1)*uj[3]+(x0*y1*z0)*uj[4]+(x1*y1*z0)*uj[5]+(x1*y1*z1)*uj[6]+(x0*y1*z1)*uj[7]; // perform trilinear interpolation
			c2 = shading(rainbow_color(255.0f*def_scale_u*length(ui)), p+p2, normal, camera_cache); // rainbow_color(255.0f*def_scale_u*length(ui));
		}
		draw_triangle_interpolated(p+p0, p+p1, p+p2, c0, c1, c2, camera_cache, bitmap, zbuffer); // draw triangle with interpolated colors
	}
}

#endif // Graphics