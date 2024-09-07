//! This module contains code related to ligation. For example, using restriction enzymes to
//! combine or otherwise edit DNA segments.

use crate::{sequence::Nucleotide, Seq};

/// Digest the sequence with one or more REs.
pub fn digest(selected: &[String], seq: &[Nucleotide]) -> Vec<Seq> {
    let mut result = Vec::new();

    result
}
