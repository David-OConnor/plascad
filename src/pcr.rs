//! This module assists in identifying PCR parameters

use bincode::{Decode, Encode};

use crate::{
    gui::navigation::{Page, PageSeq},
    primer::Primer,
    sequence::{Nucleotide, Seq},
    util::RangeIncl,
    PcrUi, State,
};

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
    pub fn extension_time(&self, product_len: usize) -> u16 {
        match self {
            Self::NormalFidelity => (60 * product_len / 1_000) as u16,
            // 15 - 30. 15 recommended in FastCloning guide.
            Self::HighFidelity => (15 * product_len / 1_000) as u16,
        }
    }

    pub fn to_str(self) -> String {
        match self {
            Self::NormalFidelity => "Normal fidelity",
            Self::HighFidelity => "High fidelity (eg Phusion)",
        }
        .to_owned()
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
            // Alternative: Ta = 0.3 x  Tm(primer) + 0.7 Tm(product) – 14.9.
            annealing: TempTime::new(data.primer_tm - 5., 30), // 15-60s. How do we choose.
            // 72 is good if Taq, and Phusion.
            extension: TempTime::new(72., data.polymerase_type.extension_time(data.product_len)),
            // Alternatively: 5-10 mins? Perhaps 30s per 1kb?)
            final_extension: TempTime::new(72., 60),
            num_cycles: data.num_cycles,
        }
    }
}

/// Create a new tab containing of the PCR amplicon.
pub fn make_amplicon_tab(
    state: &mut State,
    product_seq: Seq,
    range: RangeIncl,
    fwd_primer: Primer,
    rev_primer: Primer,
) {
    let mut product_features = Vec::new();
    for feature in &state.generic[state.active].features {
        // todo: Handle circle wraps etc.
        if range.start < feature.range.start && range.end > feature.range.end {
            let mut product_feature = feature.clone();
            // Find the indexes in the new product sequence.
            product_feature.range.start -= range.start - 1;
            product_feature.range.end -= range.start - 1;

            product_features.push(product_feature);
        }
    }

    state.add_tab();

    // if let Some(seq) = range_combined.index_seq(&state.generic.seq) {
    state.generic[state.active].seq = product_seq;

    // Include the primers used for PCR, and features that are included in the new segment.
    // note that the feature indexes will need to change.

    let product_primers = vec![fwd_primer, rev_primer];

    state.generic[state.active].features = product_features;
    state.generic[state.active].primers = product_primers;
    state.generic[state.active].metadata.plasmid_name = "PCR amplicon".to_owned();

    state.sync_seq_related(None);
    // state.sync_primer_metrics();

    state.ui.page = Page::Sequence;
    state.ui.page_seq = PageSeq::View;
}
