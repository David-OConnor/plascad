//! This module contains code related to ligation. For example, using restriction enzymes to
//! combine or otherwise edit DNA segments.

use crate::{sequence::Nucleotide, Seq};
use crate::restriction_enzyme::{ReMatch, RestrictionEnzyme};
use crate::sequence::SeqTopology;


pub struct LigationFragment {
    pub seq: Seq,
    pub re_left: RestrictionEnzyme,
    pub re_right: RestrictionEnzyme,
}

/// Digest the sequence with one or more REs.
/// `matches` here is all matches; we filter by selected here.
pub fn digest(selected: &[String], matches: &[ReMatch], re_lib: &[RestrictionEnzyme], seq: &[Nucleotide],
              topology: SeqTopology) -> Vec<LigationFragment> {
    let mut result = Vec::new();

    let mut cut_locs = Vec::new();
    for re_match in matches {
        if re_match.lib_index +1 > re_lib.len() {
            eprintln!("Error with RE lib");
            continue;
        }
        let re = &re_lib[re_match.lib_index];

        if !selected.contains(&re.name) {
            continue
        }

        cut_locs.push(re_match.seq_index);
    }

    if cut_locs.is_empty() {
        return result;
    }

    let mut cut_loc = cut_locs[0];
    let mut cuts_i = 0;
    let mut current_fragment = Vec::new();

    for seq_i in cut_loc..seq.len() {
        let nt = seq[seq_i];

        if seq_i == cut_loc {
            if !current_fragment.is_empty() {
                result.push(current_fragment.clone());
            }

            current_fragment = Vec::new();

            if cuts_i + 1 < cut_locs.len() {
                cuts_i += 1;
            }
            cut_loc = cut_locs[cuts_i];
        }
        current_fragment.push(nt);
    }

    match topology {
        SeqTopology::Circular => {
            // Create the final fragment, between the last and first cut sites.
            for seq_i in 0..cut_locs[0] {
                let nt = seq[seq_i];
                current_fragment.push(nt);

                if seq_i == cut_locs[0] {
                    break;
                    // result.push(current_fragment.clone());
                    //
                    // current_fragment = Vec::new();
                    //
                    // if cuts_i + 1 < cut_locs.len() {
                    //     cuts_i += 1;
                    // }
                    // cut_loc = cut_locs[cuts_i];
                }
            }
            // todo: This isn't right. Continue to the first cut site.
            // From the last cut site to the first, wrapping through the origin.
            result.push(current_fragment);
        }
        SeqTopology::Linear => {
            // From the origin to the first cut site.
            result.push(seq[..cut_locs[0]].to_vec());
            // From the last cut site to the end.
            result.push(current_fragment);
        }
    }

    result
}
