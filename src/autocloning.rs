//! Used for a semi-automated cloning process that chooses a suitable backbone and restriction enzymes
//! or primers.

use crate::{backbones::Backbone, restriction_enzyme::RestrictionEnzyme, sequence::Nucleotide};

/// For a given insert and vector, find suitable  restriction enzymes for cloning
/// todo: Fow now, we limit our options to unique-match, single-cutters.
/// todo: Combine this filter code with that in ligation.
pub fn find_re_candidates(
    backbone: &Backbone,
    seq_insert: &[Nucleotide],
    lib: &[RestrictionEnzyme],
) -> Vec<RestrictionEnzyme> {
    let mut result = Vec::new();

    result
}
