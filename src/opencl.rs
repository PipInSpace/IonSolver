// OpenCL functions
use ocl::{flags, Buffer, Device, Platform, Queue};

pub fn device_selection(domains: u32) -> Vec<Device> {
    let devices = get_devices();
    if devices.is_empty() {
        println!("No OpenCL Device detected. Aborting...");
        std::process::exit(1);
    } else if devices.len() == 1 {
        println!("1 OpenCL device detected");
    } else {
        println!("{} OpenCL devices detected", devices.len());
    }
    for dev in &devices {
        println!("  - {}", dev.name().expect("Device should have name"))
    }
    let mut device_infos: Vec<Device> = vec![devices[0]; domains as usize]; // Is completely overwritten
    let mut device_type_ids: Vec<Vec<Device>> = vec![]; // Device auto-selection
    for d in devices {
        let name_1 = d
            .name()
            .expect("Unable to get device name in auto-selection");
        let mut already_exists = false;
        for d_comp in &mut device_type_ids {
            let name_2 = d_comp[0].name().unwrap();
            if name_1 == name_2 {
                d_comp.push(d);
                already_exists = true;
            }
        }
        if !already_exists {
            device_type_ids.push(vec![d]);
        }
    }

    let mut best_value = 0;
    let mut best_j = -1i32;
    for j in 0..device_type_ids.len() {
        let value = get_tflops(device_type_ids[j][0]); //TODO: Better algorithm
        if device_type_ids.len() >= domains.try_into().unwrap() && value > best_value {
            best_value = value;
            best_j = j as i32;
        }
    }
    if best_j >= 0 {
        device_infos[..(domains as usize)]
            .copy_from_slice(&device_type_ids[best_j as usize][..(domains as usize)]);
        //for d in 0..domains as usize {
        //    device_infos[d] = device_type_ids[best_j as usize][d];
        //}
    } else {
        println!("Warning! Not enough devices of the same type available. Using single fastest device for all domains.");
        for d in device_infos.iter_mut().take(domains as usize) {
            *d = get_device_with_most_flops();
        }
    }
    device_infos
}

pub fn get_devices() -> Vec<Device> {
    let platform = Platform::default();

    Device::list_all(platform).expect("Cannot find devices")
}

/// Combines the graphics and simulation source files and
/// removes embedded default defines needed for syntax highlighting
pub fn get_opencl_code() -> String {
    let sim_source: Vec<&str> = include_str!("sim_kernels.cl")
        .split("EndTempDefines%")
        .collect();
    let graphics_source: Vec<&str> = include_str!("graphics_kernels.cl")
        .split("EndTempDefines%")
        .collect();
    sim_source[1].to_string() + graphics_source[1]
}

pub fn create_buffer<T: ocl::OclPrm, I: Into<ocl::SpatialDims> + Clone>(
    queue: &Queue,
    size: I,
    fill_value: T,
) -> Buffer<T> {
    if (size.clone().into() as ocl::SpatialDims).to_len() >= 1 {
        return Buffer::<T>::builder()
            .queue(queue.clone())
            .len(size)
            .fill_val(fill_value)
            .flags(flags::MEM_READ_WRITE)
            .build()
            .unwrap();
    }
    // use size of 1 if invalid
    return Buffer::<T>::builder()
        .queue(queue.clone())
        .len([1])
        .fill_val(fill_value)
        .flags(flags::MEM_READ_WRITE)
        .build()
        .unwrap();
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
    for (i, &device) in devices.iter().enumerate() {
        if get_tflops(devices[i]) > best_value {
            best_value = get_tflops(device);
            best_i = i as i32;
        }
    }
    devices[best_i as usize]
}

#[macro_export]
/// Builds a kernel with named arguments from a program, queue, kernel name and work size. Adds named arguments given as tuples of ("name", arg).
macro_rules! kernel_n {
    ($p:expr, $q:expr, $name:expr, $n:expr, $( $arg:expr),*) => {
        {
            let mut kernel_builder = Kernel::builder();
            kernel_builder.program(&$p).name($name).queue($q.clone()).global_work_size($n);
            $(
                kernel_builder.arg_named($arg.0, $arg.1);
            )*
            kernel_builder.build().unwrap()
        }
    };
}

#[macro_export]
/// Builds a kernel with unnamed arguments from a program, queue, kernel name and work size. Adds unnamed arguments.
macro_rules! kernel {
    ($p:expr, $q:expr, $name:expr, $n:expr, $( $arg:expr),*) => {
        {
            let mut kernel_builder = Kernel::builder();
            kernel_builder.program(&$p).name($name).queue($q.clone()).global_work_size($n);
            $(
                kernel_builder.arg($arg);
            )*
            kernel_builder.build().unwrap()
        }
    };
}

#[macro_export]
/// Creates a KernelBuilder from a program, queue, kernel name and work size.
macro_rules! kernel_builder {
    ($p:expr, $q:expr, $name:expr, $n:expr) => {
        {
            let mut kernel_builder = Kernel::builder();
            kernel_builder.program(&$p).name($name).queue($q.clone()).global_work_size($n);
            kernel_builder
        }
    };
}

#[macro_export]
/// Appends named arguments given as tuples of ("name", arg).
macro_rules! kernel_args_n {
    ($kernel:expr, $( $arg:expr),*) => {
        $(
            $kernel.arg_named($arg.0, $arg.1);
        )*
    };
}

#[macro_export]
/// Appends unnamed arguments.
macro_rules! kernel_args {
    ($kernel:expr, $( $arg:expr),*) => {
        $(
            $kernel.arg($arg);
        )*
    };
}