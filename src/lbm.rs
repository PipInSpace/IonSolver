//Defines

#[allow(dead_code)]
pub enum VelocitySet {
    D2Q9,
    D3Q15,
    D3Q19,
    D3Q27,
}

#[allow(dead_code)]
pub enum Extensions {

}

#[derive(Clone, Copy)]
pub struct Lbm {
    pub n_x: u32, //Size
    pub n_y: u32,
    pub n_z: u32,

    pub d_x: u32, //Domain
    pub d_y: u32,
    pub d_z: u32,

    pub o_x: u32, //Offset
    pub o_y: u32,
    pub o_z: u32,

    pub nu: f32,

    pub dimensions: u8,
    pub velocity_set: u8,
    pub transfers: u8,
}

impl Lbm {
    pub fn new(vel_set_cfg: VelocitySet) -> Lbm {
        match vel_set_cfg {
            VelocitySet::D2Q9 => Lbm {
                n_x: 1,
                n_y: 1,
                n_z: 1,
                d_x: 1,
                d_y: 1,
                d_z: 1,
                o_x: 0,
                o_y: 0,
                o_z: 0,
                nu: 1.0f32 / 6.0f32,
                dimensions: 2,
                velocity_set: 9,
                transfers: 3,
            },
            VelocitySet::D3Q15 => Lbm {
                n_x: 1,
                n_y: 1,
                n_z: 1,
                d_x: 1,
                d_y: 1,
                d_z: 1,
                o_x: 0,
                o_y: 0,
                o_z: 0,
                nu: 1.0f32 / 6.0f32,
                dimensions: 3,
                velocity_set: 15,
                transfers: 5,
            },
            VelocitySet::D3Q19 => Lbm {
                n_x: 1,
                n_y: 1,
                n_z: 1,
                d_x: 1,
                d_y: 1,
                d_z: 1,
                o_x: 0,
                o_y: 0,
                o_z: 0,
                nu: 1.0f32 / 6.0f32,
                dimensions: 3,
                velocity_set: 19,
                transfers: 5,
            },
            VelocitySet::D3Q27 => Lbm {
                n_x: 1,
                n_y: 1,
                n_z: 1,
                d_x: 1,
                d_y: 1,
                d_z: 1,
                o_x: 0,
                o_y: 0,
                o_z: 0,
                nu: 1.0f32 / 6.0f32,
                dimensions: 3,
                velocity_set: 27,
                transfers: 9,
            },
        }
    }

    pub fn get_opencl_code(self) -> String {
        return include_str!("kernels.cl").to_owned() + &self.get_device_defines();
    }

    fn get_device_defines(self) -> String {
        return
         "\n    #define def_Nx ".to_owned() + &self.n_x.to_string()+"u"
        +"\n	#define def_Ny " + &self.n_y.to_string()+"u"
        +"\n	#define def_Nz " + &self.n_z.to_string()+"u"
        +"\n	#define def_N  " + &self.get_n().to_string()+"ul"
    
        +"\n	#define def_Dx " + &self.d_x.to_string()+"u"
        +"\n	#define def_Dy " + &self.d_y.to_string()+"u"
        +"\n	#define def_Dz " + &self.d_z.to_string()+"u"
    
        +"\n	#define def_Ox " + &self.o_x.to_string()+"" // offsets are signed integer!
        +"\n	#define def_Oy " + &self.o_y.to_string()+""
        +"\n	#define def_Oz " + &self.o_z.to_string()+""
    
        +"\n	#define def_Ax "+ &(self.n_y * self.n_z).to_string()+"u"
        +"\n	#define def_Ay "+ &(self.n_z * self.n_x).to_string()+"u"
        +"\n	#define def_Az "+ &(self.n_x * self.n_y).to_string()+"u"
    
        +"\n	#define D"+ &self.dimensions.to_string()+"Q"+ &self.velocity_set.to_string()+"" // D2Q9/D3Q15/D3Q19/D3Q27
        +"\n	#define def_velocity_set "+ &self.velocity_set.to_string()+"u" // LBM velocity set (D2Q9/D3Q15/D3Q19/D3Q27)
        +"\n	#define def_dimensions "+ &self.dimensions.to_string()+"u" // number spatial dimensions (2D or 3D)
        +"\n	#define def_transfers "+ &self.transfers.to_string()+"u" // number of DDFs that are transferred between multiple domains
    
        +"\n	#define def_c 0.57735027f" // lattice speed of sound c = 1/sqrt(3)*dt
        +"\n	#define def_w " + &self.nu.to_string()+"f" // relaxation rate w = dt/tau = dt/(nu/c^2+dt/2) = 1/(3*nu+1/2)
    ;
    }

    fn get_n(self) -> u64 {
        self.n_x as u64 * self.n_y as u64 * self.n_z as u64
    }
}

fn get_tau() -> f32 {
    3.0f32 * get_nu() + 0.5f32
}
fn get_nu() -> f32 {
    1.0f32 / 6.0f32
}