//! This module assists in identifying PCR parameters

use bincode::{Decode, Encode};
use crate::PcrUi;

/// This is a common pattern for PCR parameters
#[derive(Default, Encode, Decode)]
pub struct TempTime {
    /// In °C
    pub temp: f32,
    /// In seconds
    pub time: u16,
}

impl TempTime {
    pub fn new(temp: f32, time: u16) -> Self {
        Self { temp, time }
    }
}

#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum PolymeraseType {
    NormalFidelity,
    /// Eg Phusion; results in a shorter extension time.
    HighFidelity,
}

impl Default for PolymeraseType {
    fn default() -> Self {
        Self::NormalFidelity
    }
}

impl PolymeraseType {
    pub fn extension_time(&self, product_len: usize) -> u16{
        match self {
            Self::NormalFidelity => (60 * product_len / 1_000) as u16,
            // 15 - 30. 15 recommended in FastCloning guide.
            Self::HighFidelity => (15 * product_len / 1_000) as u16,
        }
    }

    pub fn to_str(&self) -> String {
        match self {
            Self::NormalFidelity => "Normal fidelity",
            Self::HighFidelity => "High fidelity (eg Phusion)"
        }.to_owned()
    }
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

impl PcrParams {
    pub fn new(data: &PcrUi) -> Self {
        Self {
            // 94-98? 30-120s?
            initial_denaturation: TempTime::new(94., 120),
            // 94-98? 10-30s?
            denaturation: TempTime::new(94., 30),
            annealing: TempTime::new(data.primer_tm - 5., 30), // 15-60s. How do we choose.
            // 72 is good if Taq, and Phusion.
            extension: TempTime::new(72., data.polymerase_type.extension_time(data.product_len)),
            // Alternatively: 5-10 mins? Perhaps 30s per 1kb?)
            final_extension: TempTime::new(72., 60),
            num_cycles: data.num_cycles,
        }
    }
}
