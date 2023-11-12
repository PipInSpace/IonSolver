#![allow(non_snake_case)]
extern crate ocl;
use std::time::{Duration, Instant};
use std::{f32::consts::PI, fs, sync::mpsc, thread};

mod efield_precompute;
mod graphics;
mod info;
mod lbm;
mod opencl;
mod setup;
mod units;
use eframe::*;
use egui::{Color32, ColorImage, InnerResponse, Label, Sense, Stroke, TextureOptions, Vec2};
use lbm::*;

/// Is send by the simulation to communicate information to the UI
pub struct SimState {
    step: u32,
    img: ColorImage,
}

impl SimState {
    pub fn new() -> SimState {
        Self {
            step: 0,
            img: ColorImage::default(),
        }
    }
}

impl Default for SimState {
    fn default() -> Self {
        Self::new()
    }
}

/// Is send by the UI to communicate information to the simulation
pub struct SimControlTx {
    //SimControlTx is used to send information to the simulation thread
    paused: bool,
    save: bool,
    clear_images: bool,
    frame_spacing: u32,
    active: bool,
    camera_rotation: Vec<f32>,
    camera_zoom: f32,
}

struct SimControl {
    //SimControl handles control/displaying of the Simulation over the UI
    paused: bool,
    save: bool,
    step: u32,
    clear_images: bool,
    frame_spacing: u32,
    frame_spacing_str: String,
    display_img: Option<egui::TextureHandle>,
    camera_rotation: Vec<f32>,
    camera_zoom: f32,

    ctrl_tx: mpsc::Sender<SimControlTx>,
    sim_rx: mpsc::Receiver<SimState>,
}

impl SimControl {
    pub fn reset_sim(&mut self) {
        self.paused = true;
    }

    pub fn send_control(&self) {
        self.ctrl_tx
            .send(SimControlTx {
                paused: self.paused,
                save: self.save,
                clear_images: self.clear_images,
                frame_spacing: self.frame_spacing,
                active: true,
                camera_rotation: self.camera_rotation.clone(),
                camera_zoom: self.camera_zoom,
            })
            .expect("GUI cannot communicate with sim thread");
    }

    // Egui response ro
    fn draw_top_controls(
        &mut self,
        ctx: &egui::Context,
        _frame: &mut eframe::Frame,
        send_control: &mut bool,
    ) -> InnerResponse<InnerResponse<()>> {
        //egui styles
        let zeromargin = egui::Margin {
            left: 0.0,
            right: 0.0,
            top: 0.0,
            bottom: 0.0,
        };
        let small_left_margin = egui::Margin {
            left: 0.0,
            right: 0.0,
            top: 1.0,
            bottom: -1.0,
        };
        let frame = egui::Frame {
            fill: Color32::WHITE,
            inner_margin: small_left_margin,
            outer_margin: zeromargin,
            ..Default::default()
        };
        let transparent_stroke = Stroke {
            width: 0.0,
            color: Color32::TRANSPARENT,
        };

        egui::TopBottomPanel::top("top_controls")
            .frame(frame)
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    let spacing = egui::style::Spacing {
                        item_spacing: egui::Vec2 { x: 0., y: 0. },
                        ..Default::default()
                    };
                    ui.style_mut().spacing = spacing;
                    ui.style_mut().visuals.widgets.hovered.expansion = 0.0;
                    ui.style_mut().visuals.widgets.hovered.weak_bg_fill =
                        Color32::from_rgb(0xCC, 0xCC, 0xCC);
                    ui.style_mut().visuals.widgets.hovered.rounding = 0.0.into();
                    ui.style_mut().visuals.widgets.inactive.weak_bg_fill = Color32::WHITE;
                    ui.style_mut().visuals.widgets.inactive.rounding = 0.0.into();
                    ui.style_mut().visuals.widgets.active.rounding = 0.0.into();
                    ui.style_mut().visuals.widgets.noninteractive.rounding = 0.0.into();
                    if ui
                        .add_sized(
                            [40., 18.],
                            egui::Button::new(if self.paused { "Play" } else { "Pause" })
                                .rounding(0.0f32)
                                .stroke(transparent_stroke),
                        )
                        .clicked()
                    {
                        self.paused = !self.paused;
                        *send_control = true;
                    }
                    if ui
                        .add(
                            egui::Button::new("Reset")
                                .rounding(0.0f32)
                                .stroke(transparent_stroke),
                        )
                        .clicked()
                    {
                        self.reset_sim();
                    }
                    if ui
                        .add(
                            egui::Button::new(if self.save {
                                "Saving Enabled"
                            } else {
                                "Saving Disabled"
                            })
                            .rounding(0.0f32)
                            .stroke(transparent_stroke),
                        )
                        .clicked()
                    {
                        self.save = !self.save;
                        *send_control = true;
                    }
                    if ui
                        .add(
                            egui::Button::new(if self.clear_images {
                                "Old Output Deleted"
                            } else {
                                "Old Output Kept"
                            })
                            .rounding(0.0f32)
                            .stroke(transparent_stroke),
                        )
                        .clicked()
                    {
                        self.clear_images = !self.clear_images;
                        *send_control = true;
                    }
                    let mut text = self.frame_spacing_str.clone();
                    if ui
                        .add(egui::TextEdit::singleline(&mut text).desired_width(50.0))
                        .changed()
                    {
                        self.frame_spacing_str = text.clone(); // clone jik text is cbv here
                        let result = str::parse::<u32>(&text);
                        if let Ok(value) = result {
                            if value > 0 {
                                self.frame_spacing = value;
                                *send_control = true;
                            }
                        }
                    }
                })
            })
    }
}

impl App for SimControl {
    /// UI update loop
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let recieve_result = self.sim_rx.try_iter().last();
        if let Some(recieve) = recieve_result {
            if self.display_img.is_none() {
                self.display_img =
                    Some(ctx.load_texture("sim", recieve.img.clone(), Default::default()))
            }
            self.display_img
                .as_mut()
                .expect("Isn't TextureHandle")
                .set(recieve.img, TextureOptions::default());
            self.step = recieve.step;
        }

        let mut send_control = false; // determines if ui needs to send control

        //egui styles
        let zeromargin = egui::Margin {
            left: 0.0,
            right: 0.0,
            top: 0.0,
            bottom: 0.0,
        };
        let small_left_margin = egui::Margin {
            left: 0.0,
            right: 0.0,
            top: 1.0,
            bottom: -1.0,
        };
        let mut frame = egui::Frame {
            fill: Color32::WHITE,
            inner_margin: small_left_margin,
            outer_margin: zeromargin,
            ..Default::default()
        };

        //Update gui
        self.draw_top_controls(ctx, _frame, &mut send_control);

        frame.outer_margin = zeromargin;
        frame.fill = Color32::BLACK;
        egui::CentralPanel::default().frame(frame).show(ctx, |ui| {
            ui.style_mut().visuals.override_text_color = Some(Color32::WHITE);
            match &self.display_img {
                Some(img) => {
                    let response = ui.image(img.id(), ui.available_size());
                    let id = ui.id();
                    let response = ui.interact(response.rect, id, Sense::click_and_drag());
                    if response.dragged() {
                        let delta = response.drag_delta();
                        if delta != Vec2::ZERO {
                            //TODO: send data to simulation thread for graphics
                            self.camera_rotation[0] += delta.x;
                            self.camera_rotation[1] =
                                (self.camera_rotation[1] - delta.y).clamp(-90.0, 90.0);
                            send_control = true;
                        }
                    }
                }
                None => {
                    ui.style_mut()
                        .text_styles
                        .get_mut(&egui::TextStyle::Body)
                        .unwrap()
                        .size = 18.0;
                    ui.with_layout(
                        egui::Layout::centered_and_justified(egui::Direction::TopDown),
                        |ui| {
                            ui.label("Simulation Graphic Output");
                            ui.add(Label::new("text"));
                        },
                    );
                }
            }
            ui.input(|i| {
                if i.scroll_delta != Vec2::ZERO {
                    self.camera_zoom = (self.camera_zoom
                        + (self.camera_zoom * (i.scroll_delta.y * 0.001)))
                        .max(0.01);
                    send_control = true;
                }
            });
        });

        if send_control {
            self.send_control()
        }

        ctx.request_repaint_after(Duration::from_millis(100));

        thread::spawn(move || {});
    }
}

pub fn load_icon(icon_bytes: &[u8]) -> Option<eframe::IconData> {
    if let Ok(image) = image::load_from_memory(icon_bytes) {
        let image = image.to_rgba8();
        let (width, height) = image.dimensions();
        Some(eframe::IconData {
            width,
            height,
            rgba: image.as_raw().to_vec(),
        })
    } else {
        None
    }
}

fn main() {
    println!("{}", info::LOGO_COLOR);

    // UI params
    let icon_bytes = include_bytes!("../icons/IonSolver.png");
    let options = eframe::NativeOptions {
        initial_window_size: Some(egui::vec2(1000.0, 650.0)),
        follow_system_theme: false,
        default_theme: Theme::Light,
        icon_data: load_icon(icon_bytes.as_ref()),
        ..Default::default()
    };

    let (sim_tx, sim_rx) = mpsc::channel();
    let (ctrl_tx, ctrl_rx) = mpsc::channel();
    let exit_ctrl_tx = ctrl_tx.clone();
    // Start Simulation loop
    let handle = thread::spawn(move || {
        //TODO: Setup mpsc messaging for data transmission
        simloop(sim_tx, ctrl_rx);
    });

    // setup simcontrol struc to pass params and channels to GUI
    let simcontrol = SimControl {
        paused: true,
        save: false,
        step: 0,
        clear_images: true,
        frame_spacing: 100,
        frame_spacing_str: "100".to_string(),
        display_img: None,
        camera_rotation: vec![0.0; 2],
        camera_zoom: 3.0,
        ctrl_tx,
        sim_rx,
    };

    // Start window loop
    eframe::run_native(
        "IonSolver",
        options,
        Box::new(|_| Box::<SimControl>::new(simcontrol)),
    )
    .expect("unable to open window");

    _ = exit_ctrl_tx.send(SimControlTx {
        paused: true,
        save: false,
        clear_images: true,
        frame_spacing: 10,
        active: false,
        camera_rotation: vec![0.0; 2],
        camera_zoom: 3.0,
    });

    handle.join().unwrap();
}

/// Runs the simulation and it's control logic in a different thread
fn simloop(sim_tx: mpsc::Sender<SimState>, ctrl_rx: mpsc::Receiver<SimControlTx>) {
    //sim_tx.send(state) sends data to the main window loop
    let mut state = SimControlTx {
        paused: true,
        save: false,
        clear_images: true,
        frame_spacing: 10,
        active: true,
        camera_rotation: vec![0.0; 2],
        camera_zoom: 3.0,
    };
    let mut i: u32 = 0;

    let mut lbm = setup::setup();
    let now = std::time::Instant::now();
    println!("Starting initialize kernel");
    lbm.initialize();
    println!("initialize kernel took: {}", now.elapsed().as_millis());

    // get initial config from ui
    let recieve_result = ctrl_rx.try_recv();
    if let Ok(recieve) = recieve_result {
        state = recieve;
    }

    // Clearing out folder if requested
    if state.save && state.clear_images {
        match fs::remove_dir_all("out") {
            Ok(_) => (),
            Err(_) => println!("Did not find out folder. Creating it."),
        }
        fs::create_dir("out").unwrap();
    }

    // Create graphics thread with evil pointer hacks
    if lbm.config.graphics_config.graphics {
        let lbm_ptr = &lbm as *const _ as usize;
        let state_ptr = &state as *const _ as usize;
        let i_ptr = &i as *const _ as usize;
        let sim_tx_g = sim_tx.clone();
        thread::spawn(move || {
            let lbm = unsafe { &*(lbm_ptr as *const Lbm) };
            let state = unsafe { &*(state_ptr as *const SimControlTx) };
            let i = unsafe { &*(i_ptr as *const u32) };
            let mut graphics_i: u32 = 0;
            let mut cached_rot: Vec<f32> = vec![0.0; 2];
            let mut cached_zoom = 0.0;

            loop {
                let now = Instant::now();
                if !state.active {
                    println!("Exiting Graphics Loop");
                    break;
                }

                // If sim not paused, only draw frame when sim time step updated
                let frame_changed = if i > &graphics_i {
                    graphics_i = *i;
                    true
                } else {
                    false
                };

                let camera_changed = cached_rot[0] != state.camera_rotation[0]
                    || cached_rot[1] != state.camera_rotation[1]
                    || cached_zoom != state.camera_zoom;
                if camera_changed {
                    //only update camera params when camera updated
                    cached_rot = state.camera_rotation.clone();
                    cached_zoom = state.camera_zoom;
                    let mut params = graphics::camera_params_rot(
                        state.camera_rotation[0] * (PI / 180.0),
                        state.camera_rotation[1] * (PI / 180.0),
                    );
                    params[0] = state.camera_zoom;
                    for d in &lbm.domains {
                        d.graphics
                            .as_ref()
                            .expect("graphics used but not enabled")
                            .camera_params
                            .write(&params)
                            .enq()
                            .unwrap();
                    }
                }
                if (camera_changed || (!state.paused && frame_changed)) && lbm.config.graphics_config.graphics {
                    // Only draws frames, never saves them
                    lbm.draw_frame(false, sim_tx_g.clone(), i);
                }
                thread::sleep(Duration::from_millis(std::cmp::max(0, 33-(now.elapsed().as_millis() as i64)) as u64));
            }
        });
    }

    let mn = (lbm.config.n_x as u64 * lbm.config.n_y as u64 * lbm.config.n_z as u64) / 1000000;
    loop {
        //This is the master loop, cannot be paused
        if !state.paused {
            let mut j = 1;
            let jnow = std::time::Instant::now();
            loop {
                //This is the loop of the simulation. Can be paused by receiving a control message
                let recieve_result = ctrl_rx.try_iter().last();
                if let Some(recieve) = recieve_result {
                    state = recieve;
                }
                if state.paused || !state.active {
                    break;
                }

                lbm.do_time_step();

                let time_per_step = (jnow.elapsed().as_micros() / j) as u32;
                if i % 20 == 0 {print!(
                    "\rStep {}, Steps/s: {}, MLUP/s: {}",
                    i,
                    1000000 / time_per_step,
                    (mn *1000000) / time_per_step as u64
                );}
                if i % state.frame_spacing == 0 && state.save && lbm.config.graphics_config.graphics
                {
                    // Saves frames if needed
                    lbm.draw_frame(true, sim_tx.clone(), &i);
                }
                i += 1;
                j += 1;

                // ##### SLOWDOWN #####
                // #####  REMOVE  #####
                // ##### SLOWDOWN #####
                //thread::sleep(Duration::from_millis(3))
                // ##### SLOWDOWN #####
                // #####  REMOVE  #####
                // ##### SLOWDOWN #####
            }
        }
        if !state.active {
            println!("\nExiting Simulation Loop");
            break;
        }
        print!("\rStep {}, Steps/s: 0, MLUP/s: 0                    ", i);
        let recieve_result = ctrl_rx.try_recv();
        if let Ok(recieve) = recieve_result {
            state = recieve;
        }
    }
}
