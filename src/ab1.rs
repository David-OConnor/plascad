//! For operations pertaining to AB1 (Applied Biosystem's sequencing) trace sequence data.

use std::collections::HashMap;

use bincode::{Decode, Encode};
use na_seq::Seq;

/// The data structure representing AB1 data.
#[derive(Clone, Debug, Default, Encode, Decode)]
pub struct SeqRecordAb1 {
    pub id: String,
    pub name: String,
    pub description: String,
    pub sequence: Seq,
    pub sequence_user: Option<Seq>,
    pub annotations: HashMap<String, String>,
    pub quality: Option<Vec<u8>>,
    pub quality_user: Option<Vec<u8>>,
    pub peak_heights: Vec<u16>,
    /// Analyzed data, for each channel.
    pub data_ch1: Vec<u16>,
    pub data_ch2: Vec<u16>,
    pub data_ch3: Vec<u16>,
    pub data_ch4: Vec<u16>,
    /// Peak locations.
    pub peak_locations: Vec<u16>,
    /// Peak locations edited by user.
    pub peak_locations_user: Option<Vec<u16>>,
}
