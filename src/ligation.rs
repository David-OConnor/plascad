//! This module contains code related to ligation. For example, using restriction enzymes to
//! combine or otherwise edit DNA segments.

use crate::{
    restriction_enzyme::{ReMatch, RestrictionEnzyme},
    sequence::{Nucleotide, SeqTopology},
    Seq,
};

pub struct LigationFragment {
    pub seq: Seq,
    /// None if the end of a linear fragment.
    pub re_left: Option<RestrictionEnzyme>,
    pub re_right: Option<RestrictionEnzyme>,
}

/// Digest the sequence with one or more REs.
/// `matches` here is all matches; we filter by selected here.
pub fn digest(
    selected: &[String],
    matches: &[ReMatch],
    re_lib: &[RestrictionEnzyme],
    seq: &[Nucleotide],
    topology: SeqTopology,
) -> Vec<LigationFragment> {
    let mut result = Vec::new();

    // Cut index, RE. If from the start or end of the seq, it's none.
    let mut cuts = Vec::new();

    for re_match in matches {
        if re_match.lib_index + 1 > re_lib.len() {
            eprintln!("Error with RE lib");
            continue;
        }
        let re = &re_lib[re_match.lib_index];

        if !selected.contains(&re.name) {
            continue;
        }

        cuts.push((re_match.seq_index, re.clone()));
    }

    if cuts.is_empty() {
        return result;
    }

    let mut cut = &cuts[0];
    let mut cuts_i = 0;
    let mut current_fragment = Vec::new();

    for seq_i in cut.0..seq.len() {
        let nt = seq[seq_i];

        if seq_i == cut.0 {
            if !current_fragment.is_empty() {
                result.push(LigationFragment {
                    seq: current_fragment.clone(),
                    re_left: Some(cuts[cuts_i - 1].1.clone()),
                    re_right: Some(cut.1.clone()),
                });
            }

            current_fragment = Vec::new();

            if cuts_i + 1 < cuts.len() {
                cuts_i += 1;
            }
            cut = &cuts[cuts_i];
        }
        current_fragment.push(nt);
    }

    match topology {
        SeqTopology::Circular => {
            // Create the final fragment, between the last and first cut sites.
            for seq_i in 0..cuts[0].0 {
                let nt = seq[seq_i];
                current_fragment.push(nt);

                if seq_i == cuts[0].0 {
                    break;
                }
            }
            // From the last cut site to the first, wrapping through the origin.
            result.push(LigationFragment {
                seq: current_fragment,
                re_left: Some(cuts[cuts.len() - 1].1.clone()),
                re_right: Some(cuts[0].1.clone()),
            });
        }
        SeqTopology::Linear => {
            // From the origin to the first cut site.
            result.push(LigationFragment {
                seq: seq[..cuts[0].0].to_vec(),
                re_left: None,
                re_right: Some(cuts[0].1.clone()),
            });

            // From the last cut site to the end.
            result.push(LigationFragment {
                seq: current_fragment,
                re_left: Some(cuts[cuts.len() - 1].1.clone()),
                re_right: None,
            });
        }
    }

    result
}

/// Ligate multiple fragments together, in each combination.
pub fn ligate(fragments: &[LigationFragment]) -> Vec<Seq> {
    let mut result = Vec::new();

    for frag in fragments {}

    result
}
