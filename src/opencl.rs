// OpenCL functions
use ocl::{Device, Platform};

pub fn device_selection(domains: u32) -> Vec<Device> {
    let devices = get_devices();
    println!("{} OpenCL device(s) detected", devices.len());
    let mut device_infos: Vec<Device> = vec![devices[0]; domains as usize]; // Is completely overwritten
                                                                            //Device auto-selection
    let mut device_type_ids: Vec<Vec<Device>> = vec![];
    for i in 0..devices.len() {
        let name_i = devices[i]
            .name()
            .expect("Unable to get device name in auto-selection");
        let mut already_exists = false;
        for j in 0..device_type_ids.len() {
            let name_j = device_type_ids[j][0].name().unwrap();
            if name_i == name_j {
                device_type_ids[j].push(devices[i]);
                already_exists = true;
            }
        }
        if !already_exists {
            device_type_ids.push(vec![devices[i]]);
        }
    }

    let mut best_value = 0;
    let mut best_j = -1i32;
    for j in 0..device_type_ids.len() {
        let value = get_tflops(device_type_ids[j][0]); //TODO: Better algorithm
        if device_type_ids.len() >= (domains as u32).try_into().unwrap() && value > best_value {
            best_value = value;
            best_j = j as i32;
        }
    }
    if best_j >= 0 {
        for d in 0..domains as usize {
            device_infos[d] = device_type_ids[best_j as usize][d];
        }
    } else {
        println!("Warning! Not enough devices of the same type available. Using single fastest device for all domains.");
        for d in 0..domains as usize {
            device_infos[d] = get_device_with_most_flops();
        }
    }
    return device_infos;
}

pub fn get_devices() -> Vec<Device> {
    let platform = Platform::default();
    let devices = Device::list_all(platform).expect("Cannot find devices");
    devices
}

pub fn get_opencl_code() -> String {
    let string: Vec<&str> = include_str!("kernels.cl").split("EndTempDefines%").collect();
    return string[1].to_string(); // Removes embedded default defines needed for syntax highlighting etc.
}

fn get_tflops(device: Device) -> i32 {
    let freq: i32 = device
        .info(ocl::enums::DeviceInfo::MaxClockFrequency)
        .unwrap()
        .to_string()
        .parse()
        .unwrap();

    let cores: i32 = device
        .info(ocl::enums::DeviceInfo::MaxComputeUnits)
        .unwrap()
        .to_string()
        .parse()
        .unwrap();
    freq * cores
}

fn get_device_with_most_flops() -> Device {
    let devices = get_devices();
    let mut best_value = 0;
    let mut best_i = -1i32;
    for i in 0..devices.len() {
        if get_tflops(devices[i]) > best_value {
            best_value = get_tflops(devices[i]);
            best_i = i as i32;
        }
    }
    devices[best_i as usize]
}

//fn get_cores(device: Device) -> u32 {
//    let compute_units: i32 = device
//        .info(ocl::enums::DeviceInfo::MaxComputeUnits)
//        .unwrap()
//        .to_string()
//        .parse()
//        .unwrap();
//    let freq: i32 = device
//        .info(ocl::enums::DeviceInfo::MaxClockFrequency)
//        .unwrap()
//        .to_string()
//        .parse()
//        .unwrap();
//    let name = device.name().unwrap().to_lowercase();
//    let nvidia_192_cores_per_cu = name.contains("pat") || freq<1000;
//    1
//}
