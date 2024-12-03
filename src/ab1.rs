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
    pub sequence_user: Option<Seq>,
    pub sequence_base: Option<Seq>,
    pub annotations: HashMap<String, String>,
    pub quality_user: Option<Vec<u8>>,
    pub quality_base: Option<Vec<u8>>,
    pub peak_heights: Vec<u16>,
    pub height_data: Vec<Vec<u16>>,
    /// Peak locations edited by user.
    pub peak_locations_user: Option<Vec<u16>>,
    /// Peak locations as called by Basecaller
    pub peak_locations_base: Option<Vec<u16>>,
}
