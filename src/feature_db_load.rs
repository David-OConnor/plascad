//! A library of known sequences we can use to automatically add features to a sequence.
//! This will eventually load data from an online source. Currently, it includes some common items,
//! and is stored in memory.

use crate::sequence::{Feature, Nucleotide};

/// todo: Abstract this out later, if it makes sense.
fn find_his_tags(seq: &[Nucleotide]) -> Vec<Feature> {
    let mut result = Vec::new();

    result
}

pub fn find_features(seq: &[Nucleotide]) -> Vec<Feature> {
    let mut result = Vec::new();

    result.append(&mut find_his_tags(seq));

    result
}
