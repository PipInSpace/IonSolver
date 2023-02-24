use std::{fs, env::args, io};

use readformat::readf;

type Vec2i = (i64, i64);

type Poly = [Vec2i; 4];

struct IonSource {
    origin: Poly,
    direction: Vec2i,
    direction_randomness: i64,
}

struct EngineData {
    walls: Vec<Poly>,
    sources: Vec<IonSource>,
}

fn main() -> Result<(), io::Error> {
    let mut args = args().skip(1);
    let data = fs::read_to_string(args.next().expect("please provide a file as an argument"))?;

    let data: EngineData = data.into();

    Ok(())
}

impl From<String> for EngineData {
    fn from(value: String) -> Self {
        let value = readf("walls: {}\nsources: {}", value.as_str()).expect("invalid data: parse step 1");
        EngineData { 
            walls: {
                let mut r = Vec::new();
                for wall in value[0].split(", ") {
                    
                }
                r
            }, 
            sources: (),
        }
    }
}
