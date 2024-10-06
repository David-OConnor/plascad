//! This module contains info related to tags, like the 6x HIS tag.

use std::ops::RangeInclusive;

use na_seq::Nucleotide;

pub struct TagMatch {
    pub lib_index: usize,
    // todo: Experimenting with ranges vice (usize, usize)
    /// 0-based indexing.
    pub seq: RangeInclusive<usize>,
}

pub struct Tag {
    pub name: String,
    // From the 5' end.
    // pub seq: Seq, // todo: You may eventually need Vec<NucleotideGeneral>.
}

impl Tag {
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_owned(),
        }
    }
}

// T7 promoter: TAATACGACTCACTATAG
// T7 term: GCTAGTTATTGCTCAGCGG
// T7 term take 2: ctagcataaccccttggggcctctaaacgggtcttgaggggttttttg

pub fn _find_tag_matches(seq: &[Nucleotide], lib: &[Tag]) -> Vec<TagMatch> {
    let mut result = Vec::new();

    result
}

/// Load a set of common Restriction enzymes. Call this at program start, to load into a state field.
pub fn _load_tag_library() -> Vec<Tag> {
    vec![
        // todo: Hisx6+x,
        // todo: T7 promoter
        // Tag::new("AatII", vec![G, A, C, G, T, C], 4),
    ]
}
