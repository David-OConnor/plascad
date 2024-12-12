//! This module assists in identifying PCR parameters

use bincode::{Decode, Encode};
use na_seq::Seq;

use crate::{
    gui::navigation::{Page, PageSeq, Tab},
    primer::{Primer, TM_TARGET},
    state::State,
    util::RangeIncl,
};

#[derive(Clone, Encode, Decode)]
/// Variables for UI fields, for determining PCR parameters.
pub struct PcrUi {
    pub primer_tm: f32,
    pub product_len: usize,
    pub polymerase_type: PolymeraseType,
    pub num_cycles: u16,
    /// Index from primer data for the load-from-primer system. For storing dropdown state.
    pub primer_selected: usize,
    /// These are for the PCR product generation
    pub primer_fwd: usize,
    pub primer_rev: usize,
}

impl Default for PcrUi {
    fn default() -> Self {
        Self {
            primer_tm: TM_TARGET,
            product_len: 1_000,
            polymerase_type: Default::default(),
            num_cycles: 30,
            primer_selected: 0,
            primer_fwd: 0,
            primer_rev: 0,
        }
    }
}

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
            // 15 - 30. 15 recommended in FastCloning guide. PHusion manual: 15-30s/kb.
            Self::HighFidelity => (15 * product_len / 1_000) as u16, //
        }
    }

    pub fn denaturation(&self) -> TempTime {
        match self {
            Self::NormalFidelity => TempTime::new(94., 30),
            Self::HighFidelity => TempTime::new(98., 10), // pHusion manual: 5-10s.
        }
    }

    pub fn denaturation_initial(&self) -> TempTime {
        match self {
            Self::NormalFidelity => TempTime::new(94., 120),
            Self::HighFidelity => TempTime::new(98., 30),
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
            initial_denaturation: data.polymerase_type.denaturation_initial(),
            denaturation: data.polymerase_type.denaturation(),
            // Alternative: Ta = 0.3 x  Tm(primer) + 0.7 Tm(product) – 14.9.
            // 15-60s. How do we choose?. Phusion manual: 10-30s.
            annealing: TempTime::new(data.primer_tm - 5., 30),
            // 72 is good if Taq, and Phusion.
            extension: TempTime::new(72., data.polymerase_type.extension_time(data.product_len)),
            // Alternatively: 5-10 mins? Perhaps 30s per 1kb?) Phusion manual: 5-10m.
            final_extension: TempTime::new(72., 60 * 6),
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
    state.tabs_open.push(Default::default());

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
