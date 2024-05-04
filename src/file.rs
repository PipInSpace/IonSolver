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

fn decode(buffer: &[u8]) -> Result<SimValues, String> {
    println!("    Decoding File");

    // Buffer position
    let mut pos: usize = 0;
    // Check header
    let mut header = "".to_owned();
    for _ in 0..15 {
        header.push(buffer[pos] as char);
        pos += 1;
    }
    if header != "IonSolver setup" {
        println!("    Invalid format!");
        return Err("Invalid format!".to_owned());
    }

    pos += 1; // increment to skip next byte so we are 4 byte aligned again

    // Values
    let n_x = to_u32(get_next_chunk(buffer, &mut pos));
    let n_y = to_u32(get_next_chunk(buffer, &mut pos));
    let n_z = to_u32(get_next_chunk(buffer, &mut pos));
    println!("    Sim size: {}, {}, {}", n_x, n_y, n_z);

    let units_m = to_f32(get_next_chunk(buffer, &mut pos));
    let units_kg = to_f32(get_next_chunk(buffer, &mut pos));
    let units_s = to_f32(get_next_chunk(buffer, &mut pos));
    let units_c = to_f32(get_next_chunk(buffer, &mut pos));

    // Walls
    println!("    Parsing Flags...");
    let mut walls: Vec<u8> = vec![];
    for _ in 0..(n_x * n_y * n_z / (4 * 8)) {
        let chunk = get_next_chunk(buffer, &mut pos);
        for byte in chunk {
            for bit in 0..8 {
                walls.push((byte >> (7 - bit)) & 1_u8);
            }
        }
    }

    // Charges
    let len = to_u32(get_next_chunk(buffer, &mut pos));
    println!("    Parsing {} Charges...", len);
    let mut charges: Vec<(u64, f32)> = Vec::with_capacity(len as usize);
    for _ in 0..len {
        let charge = to_f32(get_next_chunk(buffer, &mut pos));
        let i1 = get_next_chunk(buffer, &mut pos);
        let i2 = get_next_chunk(buffer, &mut pos);
        charges.push((to_u64(i1, i2), charge))
    }

    // Magnets
    let len = to_u32(get_next_chunk(buffer, &mut pos));
    println!("    Parsing {} Magnets...", len);
    let mut magnets: Vec<(u64, [f32; 3])> = Vec::with_capacity(len as usize);
    for _ in 0..len {
        magnets.push((
            to_u64(
                get_next_chunk(buffer, &mut pos),
                get_next_chunk(buffer, &mut pos),
            ),
            [
                to_f32(get_next_chunk(buffer, &mut pos)),
                to_f32(get_next_chunk(buffer, &mut pos)),
                to_f32(get_next_chunk(buffer, &mut pos)),
            ],
        ))
    }

    let values = SimValues {
        n_x,
        n_y,
        n_z,
        units_m,
        units_kg,
        units_s,
        units_c,
        flags: walls,
        charges,
        magnets,
    };

    println!("Completed reading file.");

    Ok(values)
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
