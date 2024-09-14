//! Used for a semi-automated cloning process that chooses a suitable backbone and restriction enzymes
//! or primers.

use crate::{backbones::Backbone, restriction_enzyme::RestrictionEnzyme, sequence::Nucleotide, StateVolatile};
use crate::ligation::{filter_multiple_seqs, filter_unique_cutters, find_common_res};
use crate::restriction_enzyme::{find_re_matches, ReMatch};
use crate::util::RangeIncl;

/// For a given insert and vector, find suitable  restriction enzymes for cloning
/// todo: Fow now, we limit our options to unique-match, single-cutters.
/// todo: Combine this filter code with that in ligation.
pub fn find_re_candidates<'a>(
    backbone: &Backbone,
    seq_insert: &[Nucleotide],
    // re_match_sets: &[&[ReMatch]],
    // insert_tab: usize,
    // insert_range: RangeIncl,
    lib: &'a [RestrictionEnzyme],
    volatile: &[StateVolatile],
// ) -> Vec<&'a RestrictionEnzyme> {
) -> Vec<RestrictionEnzyme> {
    // Note: The first part of this function is similar to how we filter REs for digest on the digest/ligation page.

    // todo: You need to make sure the insert your passing is not just the insert itself, but
    // todo includes areas around it!
    let matches_insert = find_re_matches(seq_insert, lib);
    let matches_backbone = find_re_matches(&backbone.seq, lib);

    let re_match_set = [&matches_insert, &matches_backbone];

    // Add the insert:
    // re_match_sets.push(&volatile[*active].restriction_enzyme_matches);

    // Add the backbone:
    // re_match_sets.push(&volatile[*active].restriction_enzyme_matches);

    // Set up our initial REs: Ones that match both the backbone, and insert. (anywhere, for now)
    let mut result = find_common_res(&re_match_set, lib, true);

    // For now, we are always filtering by these, as well as setting sticky ends only.
    filter_multiple_seqs(&mut result, &re_match_set, lib);
    filter_unique_cutters(&mut result, &re_match_set, lib);

    // Filter for REs that are at an appropriate place on the insert.

    result.into_iter().map(|r| r.clone()).collect()
}
