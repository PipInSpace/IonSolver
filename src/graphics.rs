use image::{ImageBuffer, Rgb};
use ocl::{flags, Buffer, Kernel, Program, Queue};

use crate::lbm::{Lbm, LbmConfig, VelocitySet};

// Each LbmDomain renders its own frame. Frames are stitched back together in the Lbm drawFrame function.
pub struct Graphics {
    kernel_clear: Kernel,
    bitmap: Buffer<i32>,
    bitmap_host: Vec<i32>,
    zbuffer: Buffer<i32>,
    zbuffer_host: Vec<i32>,
    camera_params: Buffer<f32>,

    lbm_config: LbmConfig,
    kernel_graphics_flags: Kernel,
    kernel_graphics_flags_mc: Kernel,
    kernel_graphics_field: Kernel,
    kernel_graphics_streamline: Kernel,
    kernel_graphics_q: Kernel,

    t_last_rendered_frame: u64,
}

impl Graphics {
    pub fn new(
        lbm_config: LbmConfig,
        program: &Program,
        queue: &Queue,
        flags: &Buffer<u8>,
        u: &Buffer<f32>,
    ) -> Graphics {
        let width = lbm_config.graphics_config.camera_width;
        let height = lbm_config.graphics_config.camera_height;
        let bitmap = Buffer::<i32>::builder()
            .queue(queue.clone())
            .len([width, height])
            .fill_val(0)
            .flags(flags::MEM_READ_WRITE)
            .build()
            .unwrap();
        let zbuffer = Buffer::<i32>::builder()
            .queue(queue.clone())
            .len([width, height])
            .fill_val(0)
            .flags(flags::MEM_READ_WRITE)
            .build()
            .unwrap();
        let camera_params = Buffer::<f32>::builder()
            .queue(queue.clone())
            .len(15)
            .fill_val(0.0)
            .flags(flags::MEM_READ_WRITE)
            .build()
            .unwrap();

        let kernel_clear = Kernel::builder()
            .program(&program)
            .name("graphics_clear")
            .queue(queue.clone())
            .global_work_size(bitmap.len())
            .arg_named("bitmap", &bitmap)
            .arg_named("zbuffer", &zbuffer)
            .build()
            .unwrap();

        //Basic graphics kernels:
        let kernel_graphics_flags = Kernel::builder()
            .program(&program)
            .name("graphics_flags")
            .queue(queue.clone())
            .global_work_size([lbm_config.n_x, lbm_config.n_y, lbm_config.n_z])
            .arg_named("flags", flags)
            .arg_named("camera_params", &camera_params)
            .arg_named("bitmap", &bitmap)
            .arg_named("zbuffer", &zbuffer)
            .build()
            .unwrap();
        let kernel_graphics_flags_mc = Kernel::builder()
            .program(&program)
            .name("graphics_flags_mc")
            .queue(queue.clone())
            .global_work_size([lbm_config.n_x, lbm_config.n_y, lbm_config.n_z])
            .arg_named("flags", flags)
            .arg_named("camera_params", &camera_params)
            .arg_named("bitmap", &bitmap)
            .arg_named("zbuffer", &zbuffer)
            .build()
            .unwrap();
        let kernel_graphics_field = Kernel::builder()
            .program(&program)
            .name("graphics_field")
            .queue(queue.clone())
            .global_work_size([lbm_config.n_x, lbm_config.n_y, lbm_config.n_z])
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
                .program(&program)
                .name("graphics_streamline")
                .queue(queue.clone())
                .global_work_size([
                    lbm_config.n_x / lbm_config.graphics_config.streamline_every,
                    lbm_config.n_y / lbm_config.graphics_config.streamline_every,
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
                .program(&program)
                .name("graphics_streamline")
                .queue(queue.clone())
                .global_work_size([
                    lbm_config.n_x / lbm_config.graphics_config.streamline_every,
                    lbm_config.n_y / lbm_config.graphics_config.streamline_every,
                    lbm_config.n_z / lbm_config.graphics_config.streamline_every,
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
            .program(&program)
            .name("graphics_q")
            .queue(queue.clone())
            .global_work_size([lbm_config.n_x, lbm_config.n_y, lbm_config.n_z])
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
            bitmap_host: vec![0; (width * height) as usize],
            zbuffer,
            zbuffer_host: vec![0; (width * height) as usize],
            camera_params,

            lbm_config,
            kernel_graphics_flags,
            kernel_graphics_flags_mc,
            kernel_graphics_field,
            kernel_graphics_streamline,
            kernel_graphics_q,
            t_last_rendered_frame: 0,
        }
    }

    fn enqueue_draw_frame(&mut self) {
        // camera_params = update_camera()
        // camera_params.enqueue_write
        unsafe {
            self.kernel_clear.enq().unwrap();
            //if visualisation mode
            self.kernel_graphics_q.enq().unwrap();

            self.bitmap.read(&mut self.bitmap_host).enq().unwrap();
            self.zbuffer.read(&mut self.zbuffer_host).enq().unwrap();
        }
    }
}

#[derive(Default)]
pub struct Camera {
    pub width: u32,
    pub height: u32,
}

#[derive(Clone, Copy, Default)]
pub struct GraphicsConfig {
    pub graphics: bool,
    pub background_color: u32,
    pub camera_width: u32,
    pub camera_height: u32,
    pub u_max: f32,
    pub q_min: f32,
    pub f_max: f32,
    pub streamline_every: u32,
    pub stream_line_lenght: u32,
    pub ray_transmittance: f32,
    pub ray_colour: u32,
}

impl GraphicsConfig {
    pub fn new() -> GraphicsConfig {
        GraphicsConfig {
            graphics: true,
            background_color: 0x000000,
            camera_width: 1920,
            camera_height: 1080,
            u_max: 0.25,
            q_min: 0.0001,
            f_max: 0.002,
            streamline_every: 4,
            stream_line_lenght: 128,
            ray_transmittance: 0.25,
            ray_colour: 0x005F7F,
        }
    }
}

// draw_frame function for Lbm
impl Lbm {
    pub fn draw_frame(&mut self) -> ImageBuffer<Rgb<u8>, Vec<u8>> {
        let width = self.config.graphics_config.camera_width;
        let height = self.config.graphics_config.camera_height;
        for d in 0..self.get_domain_numbers() {
            self.domains[d].graphics.enqueue_draw_frame();
        }
        for d in 0..self.get_domain_numbers() {
            self.domains[d].queue.finish().unwrap();
        }

        let mut bitmap = self.domains[0].graphics.bitmap_host.clone();
        let mut zbuffer = self.domains[0].graphics.zbuffer_host.clone();
        for d in 1..self.get_domain_numbers() {
            let bitmap_d = &self.domains[d].graphics.bitmap_host;
            let zbuffer_d = &self.domains[d].graphics.zbuffer_host;
            for i in 0..(width * height) as usize {
                let zdi = zbuffer_d[i];
                if zdi > zbuffer[i] {
                    bitmap[i] = bitmap_d[i];
                    zbuffer[i] = zbuffer_d[i];
                }
            }
        }
        let mut frame = ImageBuffer::new(width, height);
        for (x, y, pixel) in frame.enumerate_pixels_mut() {
            let color = bitmap[(y * width + x) as usize] & 0xFFFFFF;
            let r = ((color >> 16) & 0xFF) as u8;
            let g = ((color >> 8) & 0xFF) as u8;
            let b = (color & 0xFF) as u8;
            *pixel = Rgb([r, g, b]);
        }
        frame
    }
}

pub fn get_graphics_defines(graphics_config: GraphicsConfig) -> String {
    "\n	#define GRAPHICS".to_owned()
        + "\n	#define def_background_color "
        + &graphics_config.background_color.to_string()
        + ""
        + "\n	#define def_screen_width "
        + &graphics_config.camera_width.to_string()
        + "u"
        + "\n	#define def_screen_height "
        + &graphics_config.camera_height.to_string()
        + "u"
        + "\n	#define def_scale_u "
        + &(1.0f32 / (0.57735027f32 * graphics_config.u_max)).to_string()
        + "f"
        + "\n	#define def_scale_Q_min "
        + &graphics_config.q_min.to_string()
        + "f"
        + "\n	#define def_scale_F "
        + &(1.0f32 / graphics_config.f_max).to_string()
        + "f"
        + "\n	#define def_streamline_sparse "
        + &graphics_config.streamline_every.to_string()
        + "u"
        + "\n	#define def_streamline_length "
        + &graphics_config.stream_line_lenght.to_string()
        + "u"
        + "\n	#define def_n "
        + &(1.333f32).to_string()
        + "f"
        + "\n	#define COLOR_S (127<<16|127<<8|127)"
        + "\n	#define COLOR_E (  0<<16|255<<8|  0)"
        + "\n	#define COLOR_M (255<<16|  0<<8|255)"
        + "\n	#define COLOR_T (255<<16|  0<<8|  0)"
        + "\n	#define COLOR_F (  0<<16|  0<<8|255)"
        + "\n	#define COLOR_I (  0<<16|255<<8|255)"
        + "\n	#define COLOR_0 (127<<16|127<<8|127)"
        + "\n	#define COLOR_X (255<<16|127<<8|  0)"
        + "\n	#define COLOR_Y (255<<16|255<<8|  0)"
        + "\n	#define COLOR_P (255<<16|255<<8|191)"
}
