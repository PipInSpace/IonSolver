//! # types
//! 
//! Contains and organizes types used in the lbm module.
//! These types are re-exported through the lbm module.


/// Velocity discretizations in 2D and 3D.
/// 
/// - `D2Q9`:  2D
/// - `D3Q15`: 3D low precision
/// - `D3Q19`: 3D recommended
/// - `D3Q27`: 3D highest precision
/// 
#[allow(dead_code)]
#[derive(Clone, Copy, Default, serde::Serialize, serde::Deserialize)]
pub enum VelocitySet {
    #[default]
    /// 2D
    D2Q9 = 0,
    /// 3D low precision
    D3Q15 = 1,
    /// 3D recommended
    D3Q19 = 2,
    /// 3D highest precision
    D3Q27 = 3,
}

impl VelocitySet {
    pub fn get_transfers(&self) -> usize {
        match self {
            VelocitySet::D2Q9 => 3_usize,
            VelocitySet::D3Q15 => 5_usize,
            VelocitySet::D3Q19 => 5_usize,
            VelocitySet::D3Q27 => 9_usize,
        }
    }
    pub fn get_set_values(&self) -> (u8, u8, u8) {
        match self {
            //Set dimensions/velocitys/transfers from Enum
            VelocitySet::D2Q9 => (2, 9, 3),
            VelocitySet::D3Q15 => (3, 15, 5),
            VelocitySet::D3Q19 => (3, 19, 5),
            VelocitySet::D3Q27 => (3, 27, 9),
        }
    }
}

/// LBM relaxation time type.
/// 
/// - `Srt`: Single relaxation time type, more efficient
/// - `Trt`: Two-relaxation time type, more precise
/// 
#[allow(dead_code)]
#[derive(Clone, Copy, Default, serde::Serialize, serde::Deserialize)]
pub enum RelaxationTime {
    #[default]
    /// Single relaxation time, more efficient
    Srt = 0,
    /// Two-relaxation time, more precise
    Trt = 1,
}

/// Enum for different floating-point number types used in the simulation.
///
/// Types: `FP32`,`FP16S`,`FP16C`
///
/// `FP32` represents the normal floating-point number type `f32`. It takes the most memory.
///
/// `FP16S` and `FP16C` are custom floating-point number types, represented as `u16`. They take less memory.
/// `FP16S` is recommended for best precision/memory footprint.
///
/// [Learn more at this paper about custom float types.](https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats)
#[allow(dead_code)]
#[derive(Clone, Copy, Default, serde::Serialize, serde::Deserialize)]
pub enum FloatType {
    #[default]
    /// Custom float type represented as a u16, recommended
    FP16S = 0,
    /// Custom float type represented as a u16
    FP16C = 1,
    /// Default float type
    FP32 = 2,
}

impl FloatType {
    pub fn size_of(&self) -> usize {
        match self {
            FloatType::FP16S => 2_usize,
            FloatType::FP16C => 2_usize,
            FloatType::FP32 => 4_usize,
        }
    }
}

/// Enum representing a buffer of variable [`FloatType`]
///
/// [`FloatType`]: crate::lbm::FloatType
#[derive(Clone)]
pub enum VariableFloatBuffer {
    U16(ocl::Buffer<u16>), // Buffers for variable float types
    F32(ocl::Buffer<f32>),
}

/// Enum to identify a field transfered between domain boundaries
#[derive(Clone, Copy)]
pub enum TransferField {
    Fi,
    RhoUFlags,
    Qi,
}