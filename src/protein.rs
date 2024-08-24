//! Data related to proteins, eg derived from coding region features, in conjunction with
//! reading frames.

use bincode::{Decode, Encode};

use crate::{
    amino_acids::AminoAcid,
    sequence::{find_orf_matches, Feature, FeatureType, ReadingFrame},
    State,
};

#[derive(Encode, Decode)]
pub struct Protein {
    pub feature: Feature, // todo: Consider how you want to do this. Index? Only used for creation?
    pub reading_frame: ReadingFrame,
    pub seq: Vec<AminoAcid>,
}

// todo: Move this somewhere more appropriate, then store in state_volatile.
/// Find the most likely reading frame for each coding region.
pub fn sync_cr_orf_matches(state: &mut State) {
    state.volatile.cr_orf_matches = Vec::new();

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
        let mut regions = find_orf_matches(&state.generic.seq, orf);
        region_matches.append(&mut regions);
    }

    for (i, feature) in state.generic.features.iter().enumerate() {
        if feature.feature_type != FeatureType::CodingRegion {
            continue;
        }

        // Don't show binding tags.
        let label_lower = feature.label.to_lowercase().replace(" ", "");
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

        match orf_match {
            Some(m) => {
                state.volatile.cr_orf_matches.push((i, m.clone()));
            }
            None => {}
        }
    }
}
