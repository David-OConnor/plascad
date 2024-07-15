//! This module assists in identifying PCR parameters

use bincode::{Decode, Encode};

/// This is a common pattern for PCR parameters
#[derive(Default, Encode, Decode)]
pub struct TempTime {
    /// In °C
    pub temp: f32,
    /// In seconds
    pub time: u16,
}

#[derive(Default, Encode, Decode)]
pub struct PcrParams {
    pub initial_denaturation: TempTime,
    pub denaturation: TempTime,
    pub annealing: TempTime,
    pub extension: TempTime,
    pub final_extension: TempTime,
    pub num_cycles: u16,
}
