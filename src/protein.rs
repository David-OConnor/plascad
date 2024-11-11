//! Data related to proteins, eg derived from coding region features, in conjunction with
//! reading frames.

use bincode::{Decode, Encode};
use na_seq::Nucleotide;

use crate::{
    amino_acids::{AminoAcid, CodingResult},
    external_websites::PdbData,
    misc_types::{Feature, FeatureType},
    reading_frame::{find_orf_matches, ReadingFrame, ReadingFrameMatch},
    State,
};

pub const WATER_WEIGHT: f32 = 18.015; // g/mol. We subtract these when calculating a protein's weight.

// todo: Adjust A/R. Ideally, let the user customize it.
pub const HYDROPATHY_WINDOW_SIZE: usize = 9;

#[derive(Encode, Decode)]
pub struct Protein {
    pub feature: Feature, // todo: Consider how you want to do this. Index? Only used for creation?
    pub reading_frame_match: ReadingFrameMatch,
    pub aa_seq: Vec<AminoAcid>,
    pub aa_seq_precoding: Vec<AminoAcid>,
    pub aa_seq_postcoding: Vec<AminoAcid>,
    pub weight: f32,
    pub weight_with_prepost: f32,
    // window center, value.
    /// This type is for storing graphable data used by egui_plot.
    // pub hydropath_data: Vec<[f64; 2]>,
    pub hydropath_data: Vec<(usize, f32)>,
    pub pdb_data: Vec<PdbData>,
    // Note: This is more of a UI functionality; here for now.
    pub show_hydropath: bool,
}

pub fn proteins_from_seq(
    seq: &[Nucleotide],
    features: &[Feature],
    cr_orf_matches: &[(usize, ReadingFrameMatch)],
) -> Vec<Protein> {
    let mut result = Vec::new();
    for (i, feature) in features.iter().enumerate() {
        if feature.feature_type != FeatureType::CodingRegion {
            continue;
        }
        for (j, om) in cr_orf_matches {
            if *j == i {
                if let Some(seq_orf_match_dna) = om.range.index_seq(seq) {
                    // todo: DRy with feature_db_load.
                    let len = seq_orf_match_dna.len();

                    let mut aa_seq = Vec::new();
                    // We also render AAs included in the reading frame, but not specified in the coding region feature.
                    let mut aa_seq_precoding = Vec::new();
                    let mut aa_seq_postcoding = Vec::new();

                    for i_ in 0..len / 3 {
                        let i = i_ * 3; // The ORF-modified sequence index.
                        let i_actual = i + om.range.start;

                        let nts = &seq_orf_match_dna[i..i + 3];

                        // Note: We are ignoring stop codons here.
                        if let CodingResult::AminoAcid(aa) =
                            AminoAcid::from_codons(nts.try_into().unwrap())
                        {
                            if i_actual < feature.range.start {
                                aa_seq_precoding.push(aa);
                            } else if i_actual > feature.range.end {
                                aa_seq_postcoding.push(aa);
                            } else {
                                aa_seq.push(aa);
                            }
                        }
                    }

                    // todo: Def cache these weights.
                    let weight = protein_weight(&aa_seq);
                    let weight_with_prepost = weight
                        + protein_weight(&aa_seq_precoding)
                        + protein_weight(&aa_seq_postcoding);

                    let hydropath_data = hydropathy_doolittle(&aa_seq, HYDROPATHY_WINDOW_SIZE);
                    // Convert to a format accepted by our plotting lib.
                    // let mut hydropath_data = Vec::new();
                    // for (i, val) in hydropath_data_ {
                    //     hydropath_data.push([i as f64, val as f64]);
                    // }

                    result.push(Protein {
                        feature: feature.clone(),
                        reading_frame_match: om.clone(),
                        aa_seq,
                        aa_seq_precoding,
                        aa_seq_postcoding,
                        hydropath_data,
                        weight,
                        weight_with_prepost,
                        show_hydropath: true,
                        pdb_data: Vec::new(),
                    })
                }
            }
        }
    }
    result
}

/// Calculates the Kyte-Doolittle value of Hydropathy index. Uses per-AA
/// experimentally-determined hydrophobicity values, averaged over a moving window.
///
/// https://web.expasy.org/protscale/
pub fn hydropathy_doolittle(seq: &[AminoAcid], window_size: usize) -> Vec<(usize, f32)> {
    if window_size % 2 == 0 {
        eprintln!("Window size for KD must be odd");
        return Vec::new();
    }
    let mut result = Vec::new();

    let win_div_2 = window_size / 2; // Rounds down.

    if win_div_2 - 1 >= seq.len() {
        eprintln!("Error with window size for hydropathy");
        return result;
    }

    // We'll center each window on `i`.
    for i in win_div_2..seq.len() - win_div_2 - 1 {
        let aas = &seq[i - win_div_2..i + win_div_2 + 1];
        let mut val_this_window = 0.;
        for aa in aas {
            val_this_window += aa.hydropathicity();
        }
        result.push((i, val_this_window / window_size as f32));
    }

    result
}

/// The TANGO algorithm predicts beta-sheet formation propensity, which is associated with aggregation.
/// It considers factors like amino acid sequence, solvent accessibility, and secondary structure.
pub fn tango_aggregation(seq: &[AminoAcid]) {}

/// Predicts aggregation-prone regionis by calculating an aggregation score for each AA.
/// todo: Look this up.
pub fn aggrescan(seq: &[AminoAcid]) {}

// todo: Look up SLIDER aggregation tool as well.

// todo: Visualize hydrophobic and aggreg-prone regions.

// todo: Move this somewhere more appropriate, then store in state_volatile.
/// Find the most likely reading frame for each coding region.
pub fn sync_cr_orf_matches(state: &mut State) {
    state.volatile[state.active].cr_orf_matches = Vec::new();

    // These matches are for all frames.
    let mut region_matches = Vec::new();

    for orf in [
        ReadingFrame::Fwd0,
        ReadingFrame::Fwd1,
        ReadingFrame::Fwd2,
        ReadingFrame::Rev0,
        ReadingFrame::Rev1,
        ReadingFrame::Rev2,
    ] {
        let mut regions = find_orf_matches(state.get_seq(), orf);
        region_matches.append(&mut regions);
    }

    for (i, feature) in state.generic[state.active].features.iter().enumerate() {
        if feature.feature_type != FeatureType::CodingRegion {
            continue;
        }

        // Don't show binding tags.
        let label_lower = feature.label.to_lowercase().replace(' ', "");
        if label_lower.contains("xhis") || label_lower.contains("×his") {
            continue;
        }

        // Find the best reading frame match, if there is one.
        let mut orf_match = None;
        let mut smallest_match = usize::MAX;
        for rm in &region_matches {
            if !rm.range.contains(feature.range.start) || !rm.range.contains(feature.range.end) {
                continue;
            }

            // If multiple matches contain the range, choose the smallest one.
            if rm.range.len() < smallest_match {
                smallest_match = rm.range.len();
                orf_match = Some(rm.clone());
            }
        }

        if let Some(m) = orf_match {
            state.volatile[state.active]
                .cr_orf_matches
                .push((i, m.clone()));
        }
    }
}

/// calculate protein weight in kDa.
pub fn protein_weight(seq: &[AminoAcid]) -> f32 {
    if seq.is_empty() {
        return 0.; // Avoids underflow.
    }
    let mut result = 0.;
    for aa in seq {
        result += aa.weight();
    }
    (result - WATER_WEIGHT * (seq.len() - 1) as f32) / 1_000.
}

/// Convert an amino acid sequence to string, using single-letter idents.
pub fn aa_seq_to_str(seq: &[AminoAcid]) -> String {
    let mut result = String::new();

    for aa in seq {
        result.push_str(&aa.ident_single_letter());
    }

    result
}
