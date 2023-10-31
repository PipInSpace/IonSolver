//#![allow(non_snake_case)]
extern crate ocl;
use std::time::Duration;
use std::{sync::mpsc, thread, fs, io, io::Write, f32::consts::PI};

mod graphics;
mod info;
mod lbm;
mod opencl;
mod units;
mod setup;
use eframe::*;
use egui::{Color32, ColorImage, Image, Label, Sense, Stroke, TextureOptions, Vec2};
use lbm::{LbmConfig, VelocitySet};

#[allow(dead_code)]
pub struct SimState {
    step: i32,

    //Control:
    paused: bool,
    save: bool,

    img: ColorImage,
}

impl SimState {
    pub fn new() -> SimState {
        Self {
            step: 0,
            paused: true,
            save: false,
            img: ColorImage::default(),
        }
    }

    pub fn reset_sim(&mut self) {
        self.step = 0;
        self.paused = true;
    }

    pub fn update_sim(&mut self) {
        println!("Step {}", self.step);
        self.step += 1;
    }
}

impl Default for SimState {
    fn default() -> Self {
        Self::new()
    }
}

// UI+Control
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
}

impl App for SimControl {
    /// UI update loop
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let recieve_result = self.sim_rx.try_iter().last();
        if let Some(recieve) = recieve_result {
            //self.display_img = Some(ctx.load_texture("sim", ColorImage::from_rgb([1920, 1080], &recieve.img), Default::default()));
            if self.display_img.is_none() {
                self.display_img =
                    Some(ctx.load_texture("sim", recieve.img.clone(), Default::default()))
            }
            self.display_img
                .as_mut()
                .expect("Isn't TextureHandle")
                .set(recieve.img, TextureOptions::default());
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

        let transparent_stroke = Stroke {
            width: 0.0,
            color: Color32::TRANSPARENT,
        };

        //Update gui
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
                        send_control = true;
                        let _z = self.paused;
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
                        send_control = true;
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
                        send_control = true;
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
                                send_control = true;
                            }
                        }
                    }
                })
            });
        //TODO: Display the Simulation

        frame.outer_margin = zeromargin;
        frame.fill = Color32::from_rgb(0x6C, 0x6C, 0x7C);
        egui::CentralPanel::default().frame(frame).show(ctx, |ui| {
            ui.style_mut().visuals.override_text_color = Some(Color32::WHITE);
            match &self.display_img {
                Some(img) => {
                    let response = ui.add(Image::new(img.id(), ui.available_size()));
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
                    ui.add(Label::new("Simulation Graphic Output"));
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
        let sim = SimState::new();
        simloop(sim, sim_tx, ctrl_rx);
    });

    // setup simcontrol struc to pass params and channels to GUI
    let simcontrol = SimControl {
        paused: true,
        save: false,
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

fn simloop(
    sim: SimState,
    sim_tx: mpsc::Sender<SimState>,
    ctrl_rx: mpsc::Receiver<SimControlTx>,
) {
    //TODO: Simulation here
    //sim_tx.send(sim) sends data to the main window loop
    let mut state = SimControlTx {
        paused: true,
        save: false,
        clear_images: true,
        frame_spacing: 10,
        active: true,
        camera_rotation: vec![0.0; 2],
        camera_zoom: 3.0,
    };
    let mut i = 0;

    let (mut lbm, lbm_config) = setup::setup();

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

    let mut has_commenced = false; //has the simulation started
    let cached_rot: Vec<f32> = vec![0.0; 2];
    let cached_zoom = 0.0;
    loop {
        //This is the master loop, cannot be paused
        if !state.paused {
            has_commenced = true;
            loop {
                //This is the loop of the simulation. Can be paused by receiving a control message
                let recieve_result = ctrl_rx.try_iter().last();
                if let Some(recieve) = recieve_result {
                    state = recieve;
                }
                if state.paused || !state.active {
                    break;
                }

                //Simulation commences here
                lbm.do_time_step();
                if cached_rot[0] != state.camera_rotation[0]
                    || cached_rot[1] != state.camera_rotation[1]
                    || cached_zoom != state.camera_zoom
                {
                    //only update params when camera updated
                    let mut params = graphics::camera_params_rot(
                        state.camera_rotation[0] * (PI / 180.0),
                        state.camera_rotation[1] * (PI / 180.0),
                    );
                    params[0] = state.camera_zoom;
                    for d in &lbm.domains {
                        d.graphics.camera_params.write(&params).enq().unwrap();
                    }
                }

                if i % state.frame_spacing == 0 {
                    print!("\rStep {}", i);
                    io::stdout().flush().unwrap();
                    if lbm_config.graphics_config.graphics {
                        lbm.draw_frame(state.save, state.frame_spacing, sim_tx.clone(), i);
                    }
                }
                i += 1;
            }
        }
        if state.paused && state.active {
            let mut cached_rot: Vec<f32> = vec![0.0; 2];
            let mut cached_zoom = 0.0;
            loop {
                //This is the loop of the simulation if it is paused but active. useful for displaying the simulation
                let recieve_result = ctrl_rx.try_iter().last();
                if let Some(recieve) = recieve_result {
                    state = recieve;
                }
                if !state.paused || !state.active {
                    break;
                }

                if cached_rot[0] != state.camera_rotation[0]
                    || cached_rot[1] != state.camera_rotation[1]
                    || cached_zoom != state.camera_zoom
                    || !has_commenced
                {
                    //only draw when camera updated
                    has_commenced = true;
                    cached_rot = state.camera_rotation.clone();
                    cached_zoom = state.camera_zoom;
                    let mut params = graphics::camera_params_rot(
                        state.camera_rotation[0] * (PI / 180.0),
                        state.camera_rotation[1] * (PI / 180.0),
                    );
                    params[0] = state.camera_zoom;
                    for d in &lbm.domains {
                        d.graphics.camera_params.write(&params).enq().unwrap();
                    }
                    if lbm_config.graphics_config.graphics {
                        lbm.draw_frame(state.save, state.frame_spacing, sim_tx.clone(), i);
                    }
                    thread::sleep(Duration::from_millis(33)) // about 30 FPS
                }
            }
        }
        if !state.active {
            println!("\nExiting Simulation Loop");
            break;
        }
        let recieve_result = ctrl_rx.try_recv();
        if let Ok(recieve) = recieve_result {
            state = recieve;
        }
    }
}