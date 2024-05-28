use std::{fs::File, io::Write};

use ocl_macros::{bread, bwrite};

use crate::{Lbm, LbmConfig, VelocitySet, RelaxationTime, FloatType};

pub fn read<P: AsRef<std::path::Path> + std::fmt::Display>(path: P, config: &mut LbmConfig) -> Result<Lbm, String> {
    println!("\nReading simulation state from \"{}\"", path);
    let buffer: Vec<u8> = match std::fs::read(path) {
        Ok(vec) => vec,
        Err(_) => { return Err("Could not find file.".to_owned())},
    };
    decode(&buffer, config)
}

pub fn write<P: AsRef<std::path::Path> + std::fmt::Display>(lbm: &Lbm, path: P) -> Result<(), String> {
    println!("\nWriting simulation state to \"{}\"", path);
    let buffer = encode(lbm);
    match std::fs::write(path, buffer) {
        Ok(_) => { Ok(()) },
        Err(_) => { Err("writing went wrong".to_owned()) }

    }
}

impl Lbm {
    #[allow(unused)]
    pub fn write(&self) {
        let buffer = encode(self);
        let mut simfile = File::create(format!("frame{}.ion", self.get_time_step())).unwrap();
        simfile.write_all(&buffer).unwrap();
    }
}

fn decode(buffer: &[u8], config: &mut LbmConfig) -> Result<Lbm, String> {
    let mut stream = ByteStream::from_buffer(buffer);

    let mut header = String::new();
    for _ in 0..16 {
        header.push(stream.next_char());
    }

    if !header.eq("IonSolver setup\n") {
        return Err("Invalid Format!".to_owned());
    }

    config.velocity_set = match stream.next_u8() {
        0 => VelocitySet::D2Q9,
        1 => VelocitySet::D3Q15,
        2 => VelocitySet::D3Q19,
        3 => VelocitySet::D3Q27,
        _ => VelocitySet::D2Q9, // should never happen
    };

    config.relaxation_time = match stream.next_u8() {
        0 => RelaxationTime::Srt,
        1 => RelaxationTime::Trt,
        _ => RelaxationTime::Srt, // should never happen
    };

    config.float_type = match stream.next_u8() {
        0 => FloatType::FP16C,
        1 => FloatType::FP16S,
        2 => FloatType::FP32,
        _ => FloatType::FP16C, // should never happen
    };

    config.units.m = stream.next_f32();
    config.units.kg = stream.next_f32();
    config.units.s = stream.next_f32();
    config.units.a = stream.next_f32();

    config.n_x = stream.next_u32();
    config.n_y = stream.next_u32();
    config.n_z = stream.next_u32();
    
    // fixed domains are skipped for now
    stream.next_u8();

    config.d_x = stream.next_u32();
    config.d_y = stream.next_u32();
    config.d_z = stream.next_u32();

    config.nu = stream.next_f32();

    config.f_x = stream.next_f32();
    config.f_y = stream.next_f32();
    config.f_z = stream.next_f32();

    let ext = stream.next_u8();

    config.ext_equilibrium_boudaries = (ext & 0x1) != 0;
    config.ext_volume_force = (ext & 0x2) != 0;
    config.ext_force_field = (ext & 0x4) != 0;
    config.ext_magneto_hydro = (ext & 0x8) != 0;

    config.induction_range = stream.next_u8();

    let mut lbm = Lbm::new(config.to_owned());

    let d_total = lbm.config.d_x * lbm.config.d_y * lbm.config.d_z;
    let d_n = (lbm.config.n_x / lbm.config.d_x) * (lbm.config.n_y / lbm.config.d_y) * (lbm.config.n_z / lbm.config.d_z);
    
    // flags
    for d in 0..d_total {
        let mut flags: Vec<u8> = vec![];
        for _ in 0..d_n {
            flags.push(stream.next_u8());
        }

        bwrite!(lbm.domains[d as usize].flags, flags);
    }

    // densities
    for d in 0..d_total {
        let mut densities: Vec<f32> = vec![];
        for _ in 0..d_n {
            densities.push(stream.next_f32());
        }

        bwrite!(lbm.domains[d as usize].rho, densities);
    }

    // velocities
    for d in 0..d_total {
        let mut velocities: Vec<f32> = vec![];
        for _ in 0..(d_n * 3) {
            velocities.push(stream.next_f32());
        }

        bwrite!(lbm.domains[d as usize].u, velocities);
    }

    if !lbm.config.ext_magneto_hydro {
        return Ok(lbm);
    }

    // charges
    for d in 0..d_total {
        let mut charges: Vec<f32> = vec![];
        for _ in 0..d_n {
            charges.push(stream.next_f32());
        }

        bwrite!(lbm.domains[d as usize].q.as_ref().unwrap(), charges);
    }

    // static charges
    let n_charges = stream.next_u32();

    let mut static_charges: Vec<(u64, f32)> = vec![];
    for _ in 0..n_charges {
        let mut static_charge = (0u64, 0f32);
        static_charge.0 = stream.next_u64();
        static_charge.1 = stream.next_f32();
        static_charges.push(static_charge);
    }

    lbm.charges = Some(static_charges);

    // static magnets
    let n_magnets = stream.next_u32();

    let mut static_magnets: Vec<(u64, [f32; 3])> = vec![];
    for _ in 0..n_magnets {
        let mut static_magnet = (0u64, [0f32; 3]);
        static_magnet.0 = stream.next_u64();
        static_magnet.1[0] = stream.next_f32();
        static_magnet.1[1] = stream.next_f32();
        static_magnet.1[2] = stream.next_f32();
        static_magnets.push(static_magnet);
    }

    lbm.magnets = Some(static_magnets);

    if !stream.at_end() {
        return Err("Not all data could be read, file may be corrupted.".to_owned());
    }

    Ok(lbm)
}

fn encode(lbm: &Lbm) -> Vec<u8> {
    // File structure details can be found at FILE_LAYOUT.txt
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

    let domain_count = lbm.config.d_x * lbm.config.d_y * lbm.config.d_z;
    let d_n = (lbm.config.n_x * lbm.config.n_y * lbm.config.n_z) / domain_count;

    // flags
    let mut flags: Vec<u8> = vec![];
    for d_index in 0..domain_count {
        let mut flags_temp: Vec<u8> = vec![0; d_n as usize];
        bread!(lbm.domains[d_index as usize].flags, flags_temp);
        flags.extend(flags_temp.iter());
    }

    buffer.extend(flags.iter());

    // densities
    let mut densities: Vec<f32> = vec![];
    for d_index in 0..domain_count {
        let mut dens_temp: Vec<f32> = vec![0.0f32; d_n as usize];
        bread!(lbm.domains[d_index as usize].rho, dens_temp);
        densities.extend(dens_temp.iter());
    }

    for rho in densities {
        buffer.push32(rho.to_bits());
    }

    // velocities
    let mut velocities: Vec<f32> = vec![];
    for d_index in 0..domain_count {
        let mut vel_temp: Vec<f32> = vec![0.0f32; (d_n * 3) as usize];
        bread!(lbm.domains[d_index as usize].u, vel_temp);
        velocities.extend(vel_temp.iter());
    }

    for v in velocities {
        buffer.push32(v.to_bits());
    }

    if lbm.config.ext_magneto_hydro {
        // charges
        let mut charges: Vec<f32> = vec![];
        for d_index in 0..domain_count {
            let mut charge_temp: Vec<f32> = vec![0.0f32; d_n as usize];
            bread!(lbm.domains[d_index as usize].q.as_ref().unwrap(), charge_temp);
            charges.extend(charge_temp.iter());
        }

        for q in charges {
            buffer.push32(q.to_bits());
        }
    }
    
    // static charges and magnets
    
    // static charges
    if lbm.charges.is_none() {
        buffer.push32(0);
    } else {
        let charges_temp = lbm.charges.as_ref().unwrap();
        buffer.push32(charges_temp.len() as u32);
        for charge in charges_temp {
            buffer.push64(charge.0);
            buffer.push32(charge.1.to_bits()); // charge
        }
    }

    // static magnets
    if lbm.magnets.is_none() {
        buffer.push32(0);
    } else {
        let magnets_temp = lbm.magnets.as_ref().unwrap();
        buffer.push32(magnets_temp.len() as u32);
        for magnet in magnets_temp {
            buffer.push64(magnet.0);
            buffer.push32(magnet.1[0].to_bits()); // magnetization x
            buffer.push32(magnet.1[1].to_bits()); // magnetization y
            buffer.push32(magnet.1[2].to_bits()); // magnetization z
        }
    }

    buffer
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
    buffer: Vec<u8>,
    pos: usize,
}

impl ByteStream {
    pub fn from_buffer(buffer: &[u8]) -> ByteStream {
        ByteStream{buffer: buffer.to_vec(), pos: 0}
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
            value += (self.buffer[self.pos] as u32) << (i * 8);
            self.pos += 1;
        }
        value
    }
    /// gets the next u64 from stream
    pub fn next_u64(&mut self) -> u64 {
        let mut value: u64 = 0;
        for i in 0..8 {
            value += (self.buffer[self.pos] as u64) << (i * 8);
            self.pos += 1;
        }
        value
    }
    /// gets the next f32 from stream (Little endian encoding)
    pub fn next_f32(&mut self) -> f32 {
        f32::from_le_bytes(self.next_u32().to_le_bytes())
    }
    /// gets the next char from stream
    pub fn next_char(&mut self) -> char {
        self.next_u8() as char
    }

    pub fn at_end(&self) -> bool {
        self.pos == self.buffer.len()
    }
}
