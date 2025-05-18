//! For opening the browser to NCBI BLAST, PDB etc.
//!
//! PDB Search API: https://search.rcsb.org/#search-api
//! PDB Data API: https://data.rcsb.org/#data-api

use bio_apis::ncbi;

use crate::{Selection, state::State};

/// BLAST the selected Feature, primer, or selection. Prioritize the selection.
/// This function handles extracting the sequence to BLAST from possible selections.
pub fn blast(state: &State) {
    let data = &state.generic[state.active];

    let val = match state.ui.text_selection {
        Some(sel) => {
            // Don't format sel directly, as we insert the bp count downstream for use with feature selections.
            Some((
                sel.index_seq(&data.seq),
                format!("{}, {}..{}", data.metadata.plasmid_name, sel.start, sel.end),
            ))
        }
        None => match state.ui.selected_item {
            Selection::Feature(feat_i) => {
                if feat_i >= data.features.len() {
                    eprintln!("Invalid selected feature");
                    None
                } else {
                    let feature = &data.features[feat_i];
                    Some((feature.range.index_seq(&data.seq), feature.label.clone()))
                }
            }
            Selection::Primer(prim_i) => {
                if prim_i >= data.primers.len() {
                    eprintln!("Invalid selected primer");
                    None
                } else {
                    let primer = &data.primers[prim_i];
                    Some((Some(&primer.sequence[..]), primer.name.clone()))
                }
            }
            Selection::None => None,
        },
    };

    // todo: Handle reverse.

    if let Some((seq, name)) = val {
        if let Some(s) = seq {
            ncbi::open_blast(s, &name);
        }
    }
}
