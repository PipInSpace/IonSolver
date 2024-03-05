use std::{f32::consts::PI, sync::mpsc::Sender, thread};

use image::{ImageBuffer, Rgb};
use ocl::{Buffer, Kernel, Program, Queue};

use crate::*;

#[allow(unused)]
#[derive(Clone, Copy)]
/// The vector field used in visualizations like field and streamline
pub enum VecVisMode {
    U,
    E,
    B,
    EDyn,
    BDyn,
}

#[derive(Clone, Copy)]
/// Bundles arguments for graphics initialization
pub struct GraphicsConfig {
    pub graphics_active: bool, // Activate graphics engine
    pub background_color: u32,
    pub camera_width: u32,
    pub camera_height: u32,
    pub u_max: f32,
    pub q_min: f32,
    pub f_max: f32,
    pub streamline_every: u32,
    pub stream_line_lenght: u32,
    /// The vector field used in visualizations like field and streamline
    pub vec_vis_mode: VecVisMode,

    /// Visualize a vector field as streamlines
    pub streamline_mode: bool, // Active graphics modes
    /// Visualize a vector field as individual vector lines
    pub field_mode: bool,
    /// Visualize vorticity
    pub q_mode: bool,
    pub q_field_mode: bool,
    /// Visualize simulation flags as wireframe
    pub flags_mode: bool,
    /// Visualize simulation flags with marching cubes surface reconstruction
    pub flags_surface_mode: bool,
    /// Visualize coordinate system axes
    pub axes_mode: bool,
}

impl GraphicsConfig {
    pub fn new() -> GraphicsConfig {
        GraphicsConfig {
            graphics_active: true,
            background_color: 0x000000,
            camera_width: 1920,
            camera_height: 1080,
            u_max: 0.25,
            q_min: 0.0001,
            f_max: 0.002,
            streamline_every: 4,
            stream_line_lenght: 128,
            vec_vis_mode: VecVisMode::U,

            streamline_mode: false,
            field_mode: false,
            q_mode: false,
            q_field_mode: false,
            flags_mode: false,
            flags_surface_mode: false,
            axes_mode: false,
        }
    }
}

/// LbmDomain Graphics struct used to render itself to a color and z buffer.
/// 
/// Each LbmDomain renders its own frame. Different domain frames are stitched back together in the Lbm draw_frame function.
pub struct Graphics {
    kernel_clear: Kernel,
    bitmap: Buffer<i32>,
    zbuffer: Buffer<i32>,
    pub camera_params: Buffer<f32>,

    kernel_graphics_axes: Kernel,
    kernel_graphics_flags: Kernel,
    kernel_graphics_flags_mc: Kernel,
    kernel_graphics_field: Kernel,
    kernel_graphics_streamline: Kernel,
    kernel_graphics_q: Kernel,
    kernel_graphics_q_field: Kernel,

    pub streamline_mode: bool,    // Draw streamline mode
    pub field_mode: bool,         // Draw field
    pub vec_vis_mode: VecVisMode, // What Vector to visualize
    pub q_mode: bool,             // Draw q (vorticity)
    pub q_field_mode: bool,
    pub flags_mode: bool,         // Draw flags
    pub flags_surface_mode: bool, // Draw flags (surface)
    pub axes_mode: bool,         // Draw helper axes
}

impl Graphics {
    pub fn new(
        lbm_config: &LbmConfig,
        program: &Program,
        queue: &Queue,
        flags: &Buffer<u8>,
        u: &Buffer<f32>,
    ) -> Graphics {
        let width = lbm_config.graphics_config.camera_width;
        let height = lbm_config.graphics_config.camera_height;
        let n = lbm_config.n_x as u64 * lbm_config.n_y as u64 * lbm_config.n_z as u64; //TODO: use domain size
        let bitmap = opencl::create_buffer(queue, [width, height], 0i32);
        let zbuffer = opencl::create_buffer(queue, [width, height], 0i32);
        let camera_params = opencl::create_buffer(queue, 15, 0.0f32);
        camera_params.write(&new_camera_params()).enq().unwrap();

        let kernel_clear = Kernel::builder()
            .program(program)
            .name("graphics_clear")
            .queue(queue.clone())
            .global_work_size(bitmap.len())
            .arg_named("bitmap", &bitmap)
            .arg_named("zbuffer", &zbuffer)
            .build()
            .unwrap();

        //Basic graphics kernels:
        let kernel_graphics_axes = Kernel::builder()
            .program(program)
            .name("graphics_axes")
            .queue(queue.clone())
            .global_work_size(1)
            .arg_named("camera_params", &camera_params)
            .arg_named("bitmap", &bitmap)
            .arg_named("zbuffer", &zbuffer)
            .build()
            .unwrap();
        let kernel_graphics_flags = Kernel::builder()
            .program(program)
            .name("graphics_flags")
            .queue(queue.clone())
            .global_work_size([n])
            .arg_named("flags", flags)
            .arg_named("camera_params", &camera_params)
            .arg_named("bitmap", &bitmap)
            .arg_named("zbuffer", &zbuffer)
            .build()
            .unwrap();
        let kernel_graphics_flags_mc = Kernel::builder()
            .program(program)
            .name("graphics_flags_mc")
            .queue(queue.clone())
            .global_work_size([n])
            .arg_named("flags", flags)
            .arg_named("camera_params", &camera_params)
            .arg_named("bitmap", &bitmap)
            .arg_named("zbuffer", &zbuffer)
            .build()
            .unwrap();
        let kernel_graphics_field = Kernel::builder()
            .program(program)
            .name("graphics_field")
            .queue(queue.clone())
            .global_work_size([n])
            .arg_named("flags", flags)
            .arg_named("u", u)
            .arg_named("camera_params", &camera_params)
            .arg_named("bitmap", &bitmap)
            .arg_named("zbuffer", &zbuffer)
            .arg_named("slice_mode", 0)
            .arg_named("slice_x", 0)
            .arg_named("slice_y", 0)
            .arg_named("slice_z", 0)
            .build()
            .unwrap();
        let kernel_graphics_streamline = match lbm_config.velocity_set {
            VelocitySet::D2Q9 => Kernel::builder()
                .program(program)
                .name("graphics_streamline")
                .queue(queue.clone())
                .global_work_size([
                    (lbm_config.n_x / lbm_config.graphics_config.streamline_every) as u64
                        * (lbm_config.n_y / lbm_config.graphics_config.streamline_every) as u64,
                ])
                .arg_named("flags", flags)
                .arg_named("u", u)
                .arg_named("camera_params", &camera_params)
                .arg_named("bitmap", &bitmap)
                .arg_named("zbuffer", &zbuffer)
                .arg_named("slice_mode", 0)
                .arg_named("slice_x", 0)
                .arg_named("slice_y", 0)
                .arg_named("slice_z", 0)
                .build()
                .unwrap(),
            _ => Kernel::builder()
                .program(program)
                .name("graphics_streamline")
                .queue(queue.clone())
                .global_work_size([
                    (lbm_config.n_x / lbm_config.graphics_config.streamline_every) as u64
                        * (lbm_config.n_y / lbm_config.graphics_config.streamline_every) as u64
                        * (lbm_config.n_z / lbm_config.graphics_config.streamline_every) as u64,
                ])
                .arg_named("flags", flags)
                .arg_named("u", u)
                .arg_named("camera_params", &camera_params)
                .arg_named("bitmap", &bitmap)
                .arg_named("zbuffer", &zbuffer)
                .arg_named("slice_mode", 0)
                .arg_named("slice_x", 0)
                .arg_named("slice_y", 0)
                .arg_named("slice_z", 0)
                .build()
                .unwrap(),
        };
        let kernel_graphics_q = Kernel::builder()
            .program(program)
            .name("graphics_q")
            .queue(queue.clone())
            .global_work_size([n]) //TODO: this is incorrect, need own dimension size
            .arg_named("flags", flags)
            .arg_named("u", u)
            .arg_named("camera_params", &camera_params)
            .arg_named("bitmap", &bitmap)
            .arg_named("zbuffer", &zbuffer)
            .build()
            .unwrap();
        let kernel_graphics_q_field = Kernel::builder()
            .program(program)
            .name("graphics_q_field")
            .queue(queue.clone())
            .global_work_size([n]) //TODO: this is incorrect, need own dimension size
            .arg_named("flags", flags)
            .arg_named("u", u)
            .arg_named("camera_params", &camera_params)
            .arg_named("bitmap", &bitmap)
            .arg_named("zbuffer", &zbuffer)
            .build()
            .unwrap();

        Graphics {
            kernel_clear,
            bitmap,
            zbuffer,
            camera_params,

            kernel_graphics_axes,
            kernel_graphics_flags,
            kernel_graphics_flags_mc,
            kernel_graphics_field,
            kernel_graphics_streamline,
            kernel_graphics_q,
            kernel_graphics_q_field,
            vec_vis_mode: lbm_config.graphics_config.vec_vis_mode,
            streamline_mode: lbm_config.graphics_config.streamline_mode,
            field_mode: lbm_config.graphics_config.field_mode,
            q_mode: lbm_config.graphics_config.q_mode,
            q_field_mode: lbm_config.graphics_config.q_field_mode,
            flags_mode: lbm_config.graphics_config.flags_mode,
            flags_surface_mode: lbm_config.graphics_config.flags_surface_mode,
            axes_mode: lbm_config.graphics_config.axes_mode,
        }
    }
}

// draw_frame function for Lbm
impl Lbm {
    pub fn draw_frame(&self, save: bool, sim_tx: Sender<SimState>, i: &u32) {
        let width = self.config.graphics_config.camera_width;
        let height = self.config.graphics_config.camera_height;
        let domain_numbers = self.get_domain_numbers();
        let mut bitmap: Vec<i32> = vec![0; (width * height) as usize]; // Base bitmap
        let mut zbuffer: Vec<i32> = vec![0; (width * height) as usize];
        let mut bitmaps: Vec<Vec<i32>> =
            vec![vec![0; (width * height) as usize]; domain_numbers - 1]; // Holds later domain bitmaps
        let mut zbuffers: Vec<Vec<i32>> =
            vec![vec![0; (width * height) as usize]; domain_numbers - 1];
        for d in 0..domain_numbers {
            self.domains[d].enqueue_draw_frame();
        }
        for d in 0..domain_numbers {
            self.domains[d].queue.finish().unwrap();

            if d == 0 {
                self.domains[d]
                    .graphics
                    .as_ref()
                    .expect("graphics not enabled")
                    .bitmap
                    .read(&mut bitmap)
                    .enq()
                    .unwrap();
                self.domains[d]
                    .graphics
                    .as_ref()
                    .expect("graphics not enabled")
                    .zbuffer
                    .read(&mut zbuffer)
                    .enq()
                    .unwrap();
            } else {
                self.domains[d]
                    .graphics
                    .as_ref()
                    .expect("graphics not enabled")
                    .bitmap
                    .read(&mut bitmaps[d])
                    .enq()
                    .unwrap();
                self.domains[d]
                    .graphics
                    .as_ref()
                    .expect("graphics not enabled")
                    .zbuffer
                    .read(&mut zbuffers[d])
                    .enq()
                    .unwrap();
            }
        }
        self.finish_queues();

        let i = *i;
        thread::spawn(move || {
            // Generating images needs own thread for performance reasons
            for d in 0..domain_numbers - 1 {
                let bitmap_d = &bitmaps[d];
                let zbuffer_d = &zbuffers[d];
                for i in 0..(width * height) as usize {
                    let zdi = zbuffer_d[i];
                    if zdi > zbuffer[i] {
                        bitmap[i] = bitmap_d[i];
                        zbuffer[i] = zbuffer_d[i];
                    }
                }
            }

            let mut save_buffer: Vec<u8> = Vec::with_capacity(bitmap.len());

            #[cfg(not(feature = "headless"))]
            let mut pixels: Vec<egui::Color32> = Vec::with_capacity(bitmap.len());

            for pixel in &bitmap {
                let color = pixel & 0xFFFFFF;

                #[cfg(not(feature = "headless"))]
                pixels.push(egui::Color32::from_rgb(
                    ((color >> 16) & 0xFF) as u8,
                    ((color >> 8) & 0xFF) as u8,
                    (color & 0xFF) as u8,
                ));
                if save {
                    // only update save buffer if required
                    save_buffer.push(((color >> 16) & 0xFF) as u8);
                    save_buffer.push(((color >> 8) & 0xFF) as u8);
                    save_buffer.push((color & 0xFF) as u8);
                }
            }

            #[cfg(not(feature = "headless"))]
            let color_image = egui::ColorImage {
                size: [width as usize, height as usize],
                pixels,
            };

            _ = sim_tx.send(SimState {
                step: i,
                #[cfg(not(feature = "headless"))]
                img: color_image,
            }); // This may fail if the simulation is terminated, but a frame is still being generated. Error can be ignored.

            if save {
                thread::spawn(move || {
                    //Saving needs own thread for performance reasons
                    let imgbuffer: ImageBuffer<Rgb<u8>, _> =
                        ImageBuffer::from_raw(width, height, save_buffer).unwrap();
                    imgbuffer.save(format!(r"out/frame_{}.png", i)).unwrap();
                });
            }
        });
    }
}

// enqueue_draw_frame function for LbmDomain
impl LbmDomain {
    fn enqueue_draw_frame(&self) {
        let graphics = self
            .graphics
            .as_ref()
            .expect("Graphics used but not initialized");
        // Kernel enqueueing is unsafe
        unsafe {
            graphics.kernel_clear.enq().unwrap();
            if graphics.axes_mode {
                graphics.kernel_graphics_axes.enq().unwrap();
            }
            if graphics.streamline_mode {
                // Streamlines can show velocity, dynamic and static E and B field
                graphics.kernel_graphics_streamline.set_arg("u", match graphics.vec_vis_mode {
                    VecVisMode::U => &self.u,
                    VecVisMode::E => self.e.as_ref().expect("E buffer used but not initialized"),
                    VecVisMode::B => self.b.as_ref().expect("B buffer used but not initialized"),
                    VecVisMode::EDyn => self.e_dyn.as_ref().expect("E_dyn buffer used but not initialized"),
                    VecVisMode::BDyn => self.b_dyn.as_ref().expect("B_dyn buffer used but not initialized"),
                }).unwrap();
                graphics.kernel_graphics_streamline.enq().unwrap();
            }
            if graphics.field_mode {
                graphics.kernel_graphics_field.set_arg("u", match graphics.vec_vis_mode {
                    VecVisMode::U => &self.u,
                    VecVisMode::E => self.e.as_ref().expect("E buffer used but not initialized"),
                    VecVisMode::B => self.b.as_ref().expect("B buffer used but not initialized"),
                    VecVisMode::EDyn => self.e_dyn.as_ref().expect("E_dyn buffer used but not initialized"),
                    VecVisMode::BDyn => self.b_dyn.as_ref().expect("B_dyn buffer used but not initialized"),
                }).unwrap();
                graphics.kernel_graphics_field.enq().unwrap();
            }
            if graphics.q_mode {
                graphics.kernel_graphics_q.enq().unwrap();
            }
            if graphics.q_field_mode {
                graphics.kernel_graphics_q_field.enq().unwrap();
            }
            if graphics.flags_mode {
                graphics.kernel_graphics_flags.enq().unwrap();
            }
            if graphics.flags_surface_mode {
                graphics.kernel_graphics_flags_mc.enq().unwrap();
            }
        }
    }
}

#[derive(Default)]
pub struct Camera {
    pub width: u32,
    pub height: u32,
}

#[rustfmt::skip]
/// Returns a string of OpenCL C `#define`s from the provided arguments that are appended to the base OpenCl code at runtime.
pub fn get_graphics_defines(graphics_config: GraphicsConfig) -> String {
    "\n	#define GRAPHICS".to_owned()
    +"\n	#define def_background_color "  + &graphics_config.background_color.to_string()
    +"\n	#define def_screen_width "      + &graphics_config.camera_width.to_string()+"u"
    +"\n	#define def_screen_height "     + &graphics_config.camera_height.to_string()+"u"
    +"\n	#define def_scale_u "           + &(1.0f32 / (0.57735027f32 * graphics_config.u_max)).to_string()+"f"
    +"\n	#define def_scale_Q_min "       + &graphics_config.q_min.to_string()+"f"
    +"\n	#define def_scale_F "           + &(1.0f32 / graphics_config.f_max).to_string()+"f"
    +"\n	#define def_streamline_sparse " + &graphics_config.streamline_every.to_string()+"u"
    +"\n	#define def_streamline_length " + &graphics_config.stream_line_lenght.to_string()+"u"
    +"\n	#define def_n "+ &(1.333f32).to_string()+"f"
    +"\n	#define COLOR_S (127<<16|127<<8|127)"
    +"\n	#define COLOR_E (  0<<16|255<<8|  0)"
    +"\n	#define COLOR_M (255<<16|  0<<8|255)"
    +"\n	#define COLOR_T (255<<16|  0<<8|  0)"
    +"\n	#define COLOR_F (  0<<16|  0<<8|255)"
    +"\n	#define COLOR_I (  0<<16|255<<8|255)"
    +"\n	#define COLOR_0 (127<<16|127<<8|127)"
    +"\n	#define COLOR_X (255<<16|127<<8|  0)"
    +"\n	#define COLOR_Y (255<<16|255<<8|  0)"
    +"\n	#define COLOR_P (255<<16|255<<8|191)"
}

pub fn new_camera_params() -> Vec<f32> {
    let mut params: Vec<f32> = vec![0.0; 15];
    //Defaults from FluidX3D: graphics.hpp:20
    let eye_dist: u32 = 0;

    let rx = 0.5 * PI;
    let sinrx = rx.sin();
    let cosrx = rx.cos();
    let ry = PI;
    let sinry = ry.sin();
    let cosry = ry.cos();

    params[0] = 3.0; //zoom
    params[1] = 512.0; //distance from rotation center
                       //2-4 is pos x y z
                       //5-13 is a rotation matrix
    params[5] = cosrx; //Rxx
    params[6] = sinrx; //Rxy
    params[7] = 0.0; //Rxz
    params[8] = sinrx * sinry; //Ryx
    params[9] = -cosrx * sinry; //Ryy
    params[10] = cosry; //Ryz
    params[11] = -sinrx * cosry; //Rzx
    params[12] = cosrx * cosry; //Rzy
    params[13] = sinry; //Rzz
    params[14] = ((false as u32) << 31 | (false as u32) << 30 | (eye_dist & 0xFFFF)) as f32;
    params
}

pub fn camera_params_rot(rx: f32, ry: f32) -> Vec<f32> {
    let mut params: Vec<f32> = vec![0.0; 15];
    //Defaults from FluidX3D: graphics.hpp:20
    let eye_dist: u32 = 0;

    //let rx = 0.5 * PI;
    let sinrx = rx.sin();
    let cosrx = rx.cos();
    //let ry = PI;
    let sinry = ry.sin();
    let cosry = ry.cos();

    params[0] = 3.0; //zoom
    params[1] = 512.0; //distance from rotation center
                       //2-4 is pos x y z
                       //5-13 is a rotation matrix
    params[5] = cosrx; //Rxx
    params[6] = sinrx; //Rxy
    params[7] = 0.0; //Rxz
    params[8] = sinrx * sinry; //Ryx
    params[9] = -cosrx * sinry; //Ryy
    params[10] = cosry; //Ryz
    params[11] = -sinrx * cosry; //Rzx
    params[12] = cosrx * cosry; //Rzy
    params[13] = sinry; //Rzz
    params[14] = ((false as u32) << 31 | (false as u32) << 30 | (eye_dist & 0xFFFF)) as f32;
    params
}
