use std::{borrow::Borrow, fs::File, io::Write};

use crate::Lbm;

pub struct SimValues {
    pub n_x: u32,
    pub n_y: u32,
    pub n_z: u32,
    pub units_m: f32,
    pub units_kg: f32,
    pub units_s: f32,
    pub units_c: f32,
    pub flags: Vec<u8>,
    pub charges: Vec<(u64, f32)>,
    pub magnets: Vec<(u64, [f32; 3])>,
}

pub fn read<P: AsRef<std::path::Path> + std::fmt::Display>(path: P) -> Result<SimValues, String> {
    println!("Reading simulation state from \"{}\"", path);
    let buffer: Vec<u8> = std::fs::read(path).expect("Location should exist");
    decode(&buffer)
}

impl Lbm {
    #[allow(unused)]
    pub fn write(&self) {
        let buffer = encode(self);
        let mut simfile = File::create(format!("frame{}.ion", self.get_time_step())).unwrap();
        simfile.write_all(&buffer).unwrap();
    }
}

fn decode(buffer: &[u8]) -> Result<Lbm, String> {
    
}

fn encode(lbm: &Lbm) -> Vec<u8> {
    // File structure details can be found at https://github.com/PipInSpace/ionsolver-files
    let mut buffer: Vec<u8> = Vec::with_capacity(1);

    // always present and same length
    buffer.pushname(); // Header, Human-Readable
    buffer.push(lbm.config.velocity_set as u8); // Enums as u8
    buffer.push(lbm.config.relaxation_time as u8);
    buffer.push(lbm.config.float_type as u8);
    buffer.push32(lbm.config.units.m.to_bits()); // Units
    buffer.push32(lbm.config.units.kg.to_bits());
    buffer.push32(lbm.config.units.s.to_bits());
    buffer.push32(lbm.config.units.a.to_bits());
    buffer.push32(lbm.config.n_x); // Simulation sizes
    buffer.push32(lbm.config.n_y);
    buffer.push32(lbm.config.n_z);
    buffer.push(0); // Fixed domain sizes. Should the following domain sizes be used or dynamically generated?
    buffer.push32(lbm.config.d_x); // Domain sizes
    buffer.push32(lbm.config.d_y);
    buffer.push32(lbm.config.d_z);
    buffer.push32(lbm.config.nu.to_bits()); // Nu
    buffer.push32(lbm.config.f_x.to_bits()); // Volume force
    buffer.push32(lbm.config.f_y.to_bits());
    buffer.push32(lbm.config.f_z.to_bits());
    let ext: u8 = lbm.config.ext_equilibrium_boudaries as u8
        + ((lbm.config.ext_volume_force as u8) << 1)
        + ((lbm.config.ext_force_field as u8) << 2)
        + ((lbm.config.ext_magneto_hydro as u8) << 3);
    buffer.push(ext);
    buffer.push(lbm.config.induction_range);

    // values relating to simulation volume
    // domains will be saved regardless of actual coordinates to simplify
    // reading is done in the same order so it doesn't matter

    let n = lbm.config.n_x * lbm.config.n_y * lbm.config.n_z;
    let domain_count = lbm.config.d_x * lbm.config.d_y * lbm.config.d_z;
    let domain_n = lbm.domains[0].n_x * lbm.domains[0].n_y * lbm.domains[0].n_z; // since all domains are the same size this never changes

    // flags
    let mut flags: Vec<u8> = vec![];
    for d_index in 0..domain_count {
        let mut flags_temp: Vec<u8>;
        lbm.domains[d_index as usize].flags.read(&flags_temp).enq().unwrap();
        flags.extend(flags_temp.iter());
    }

    buffer.extend(flags.iter());

    // densities
    let mut densities: Vec<f32> = vec![];
    for d_index in 0..domain_count {
        let mut dens_temp: Vec<f32>;
        lbm.domains[d_index as usize].rho.read(&dens_temp).enq().unwrap();
        densities.extend(dens_temp.iter());
    }

    for rho in densities {
        buffer.push32(rho.to_bits());
    }

    // charges
    let mut charges: Vec<f32> = vec![];
    for d_index in 0..domain_count {
        let mut charge_temp: Vec<f32>;
        lbm.domains[d_index as usize].q.read(&charge_temp).enq().unwrap();
        charges.extend(charge_temp.iter());
    }

    for q in charges {
        buffer.push32(q.to_bits());
    }

    // velocities
    let mut velocities: Vec<f32> = vec![];
    for d_index in 0..domain_count {
        let mut vel_temp: Vec<f32>;
        lbm.domains[d_index as usize].u.read(&vel_temp).enq().unwrap();
        velocities.extend(vel_temp.iter());
    }

    for v in velocities {
        buffer.push32(v.to_bits());
    }

    // static charges and magnets
    
    // static charges
    if lbm.charges.is_none() {
        buffer.push32(0);
    } else {
        let charges_temp = lbm.charges.unwrap();
        buffer.push32(charges_temp.len() as u32);
        for charge in charges_temp {
            buffer.push32(charge.0 & 0x00000000FFFFFFFF); // low 4 bytes n
            buffer.push32(charge.0 >> 32); // high 4 bytes n
            buffer.push32(charge.1.to_bits()); // charge
        }
    }

    // static magnets
    if lbm.magnets.is_none() {
        buffer.push32(0);
    } else {
        let magnets_temp = lbm.magnets.unwrap();
        buffer.push32(magnets_temp.len() as u32);
        for magnet in magnets_temp {
            buffer.push32(magnet.0 & 0x00000000FFFFFFFF); // low 4 bytes n
            buffer.push32(magnet.0 >> 32); // high 4 bytes n
            buffer.push32(magnet.1.0.to_bits()); // magnetization x
            buffer.push32(magnet.1.1.to_bits()); // magnetization y
            buffer.push32(magnet.1.2.to_bits()); // magnetization z
        }
    }

    buffer
}

fn get_next_chunk(buffer: &[u8], pos: &mut usize) -> [u8; 4] {
    let mut v = [0; 4];
    v[0] = buffer[*pos];
    v[1] = buffer[*pos + 1];
    v[2] = buffer[*pos + 2];
    v[3] = buffer[*pos + 3];
    *pos += 4;

    v
}

fn to_u32(v: [u8; 4]) -> u32 {
    v[0] as u32 + ((v[1] as u32) << 8) + ((v[2] as u32) << 16) + ((v[3] as u32) << 24)
}

fn to_u64(vlow: [u8; 4], vhigh: [u8; 4]) -> u64 {
    to_u32(vlow) as u64 + ((to_u32(vhigh) as u64) << 32)
}

fn to_f32(v: [u8; 4]) -> f32 {
    f32::from_le_bytes(v)
}

/// Push different byte sizes to 8-Bit Buffer
pub trait ByteBuffer {
    fn push32(&mut self, x: u32);
    fn push64(&mut self, x: u64);
    fn pushname(&mut self);
}

impl ByteBuffer for Vec<u8> {
    /// Pushes the human-readable fileheader
    fn pushname(&mut self) {
        self.push(b'I');
        self.push(b'o');
        self.push(b'n');
        self.push(b'S');
        self.push(b'o');
        self.push(b'l');
        self.push(b'v');
        self.push(b'e');
        self.push(b'r');
        self.push(b' ');
        self.push(b's');
        self.push(b'e');
        self.push(b't');
        self.push(b'u');
        self.push(b'p');
        self.push(b'\n');
    }
    /// Push 32 bits
    fn push32(&mut self, x: u32) {
        for i in 0..4 {
            self.push(((x >> (i * 8)) & 0xFF).try_into().unwrap());
        }
    }
    /// Push 64 bits
    fn push64(&mut self, x: u64) {
        for i in 0..8 {
            self.push(((x >> (i * 8)) & 0xFF).try_into().unwrap());
        }
    }
}

pub struct ByteStream {
    buffer: &[u8],
    pos: usize,
}

impl ByteStream {
    pub fn from_buffer(buffer: &[u8]) -> ByteStream {
        ByteBuffer{buffer, 0}
    }

    /// gets the next byte from stream
    pub fn next_u8(&mut self) -> u8 {
        let value = self.buffer[self.pos];
        self.pos += 1;
        value
    }
    /// gets the next u32 from stream
    pub fn next_u32(&mut self) -> u32 {
        let mut value: u32 = 0;
        for i in 0..4 {
            value += self.buffer[self.pos] << (i * 8);
            self.pos += 1;
        }
        value
    }
    /// gets the next u64 from stream
    pub fn next_u64(&mut self) -> u64 {
        let mut value: u64 = 0;
        for i in 0..8 {
            value += self.buffer[self.pos] << (i * 8);
            self.pos += 1;
        }
        value
    }
    /// gets the next f32 from stream (Little endian encoding)
    pub fn next_f32(&mut self) -> f32 {
        f32::from_le_bytes(self.next_u32())
    }
}
