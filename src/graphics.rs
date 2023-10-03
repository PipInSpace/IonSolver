use ocl::{flags, Buffer, Kernel, Program, Queue};

use crate::lbm::{LbmConfig, VelocitySet, Lbm};

// Each LbmDomain renders its own frame. Frames are stitched back together in the Lbm drawFrame function.
pub struct Graphics {
    camera: Camera,
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
        camera: Camera,
        lbm_config: LbmConfig,
        program: &Program,
        queue: &Queue,
        flags: &Buffer<u8>,
        u: &Buffer<f32>,
        streamlines_every: u32, // Show streamlines every x cells
    ) -> Graphics {
        let bitmap = Buffer::<i32>::builder()
            .queue(queue.clone())
            .len([camera.width, camera.height])
            .fill_val(0)
            .flags(flags::MEM_READ_WRITE)
            .build()
            .unwrap();
        let zbuffer = Buffer::<i32>::builder()
            .queue(queue.clone())
            .len([camera.width, camera.height])
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
                    lbm_config.n_x / streamlines_every,
                    lbm_config.n_y / streamlines_every,
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
                    lbm_config.n_x / streamlines_every,
                    lbm_config.n_y / streamlines_every,
                    lbm_config.n_z / streamlines_every,
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
            .arg_named("camera_params", &camera_params)
            .arg_named("bitmap", &bitmap)
            .arg_named("zbuffer", &zbuffer)
            .build()
            .unwrap();

        Graphics {
            camera,
            kernel_clear,
            bitmap,
            bitmap_host: vec![],
            zbuffer,
            zbuffer_host: vec![],
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

    fn enqueue_draw_frame(&mut self){
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
    width: u32,
    height: u32,
}

// draw_frame function for Lbm
impl Lbm {
    pub fn draw_frame(&mut self) -> Vec<i32> {
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
            for i in 0..(self.domains[0].graphics.camera.width * self.domains[0].graphics.camera.height) as usize {
                let zdi = zbuffer_d[i];
                if zdi>zbuffer[i] {
                    bitmap[i] = bitmap_d[i];
                    zbuffer[i] = zbuffer_d[i];
                }
            }
        }
        return bitmap;
    }
}