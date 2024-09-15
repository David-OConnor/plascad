//! Used for a semi-automated cloning process that chooses a suitable backbone and restriction enzymes
//! or primers.

// todo: For His tags, ensure they are in frame.

use crate::{
    backbones::Backbone,
    ligation::{filter_multiple_seqs, filter_unique_cutters, find_common_res},
    restriction_enzyme::{find_re_matches, RestrictionEnzyme},
    sequence::Nucleotide,
    StateVolatile,
};

/// Include this many nucleotides to the left, and right of each insert, when searching for RE sites.
/// note: 4-6 nucleotides may be an ideal buffer. Note that this should be conservatively long.
pub const RE_INSERT_BUFFER: usize = 22;

/// For a given insert and vector, find suitable  restriction enzymes for cloning.
/// Make sure that the insert sequence is properly buffered to allow for RE matching outside
/// of the coding region (etc)'s range, upstream of this.
/// todo: Fow now, we limit our options to unique-match, single-cutters.
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

    let matches_insert = find_re_matches(seq_insert, lib);
    let matches_backbone = find_re_matches(&backbone.seq, lib);

    let re_match_set = [&matches_insert, &matches_backbone];

    // Set up our initial REs: Ones that match both the backbone, and insert. (anywhere, for now)
    let mut result = find_common_res(&re_match_set, lib, true);

    // For now, we are always filtering by these, as well as setting sticky ends only.
    filter_multiple_seqs(&mut result, &re_match_set, lib);
    filter_unique_cutters(&mut result, &re_match_set, lib);

    // Filter for REs that are at an appropriate place on the insert.

    result.into_iter().map(|r| r.clone()).collect()
}
