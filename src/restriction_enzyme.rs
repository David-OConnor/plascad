//! This module contains info related to Restriction Enzyme sites.
//!
//! [Wikipedia: List of RE sites](https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites:_A)
//! [NEB guide](https://www.neb.com/en-us/tools-and-resources/selection-charts/frequencies-of-restriction-sites)
//!
//! Note: This module only currently includes a selection of popular REs, and only ones that match
//! exact NTs.

use crate::{
    primer::PrimerDirection,
    sequence::{
        Nucleotide::{A, C, G, T},
        Seq,
    },
};

/// Unlike `Nucleotide`, this includes wildcards
#[derive(Clone, Copy)]
pub enum NucleotideGeneral {
    A,
    T,
    C,
    G,
    /// Any
    N,
    /// A or T
    W,
    /// C or G
    S,
    /// Pyrimidines: C or T
    Y,
    /// Purines: A or G
    R,
    /// A or C
    M,
    /// G or T
    K,
}

pub struct ReMatch {
    pub lib_index: usize,
    pub seq_index: usize,
    /// Direction helps align the cut site.
    pub direction: PrimerDirection,
}

pub struct RestrictionEnzyme {
    pub name: String,
    /// From the 5' end.
    pub seq: Seq, // todo: You may eventually need Vec<NucleotideGeneral>.
    /// Index to cut after, from the 5' end. For blunt ends, this will be
    /// halfway through the seq (rounded down)
    pub cut_after: u8,
}

impl RestrictionEnzyme {
    pub fn new(name: &str, seq: Seq, cut_after: u8) -> Self {
        Self {
            name: name.to_owned(),
            seq,
            cut_after,
        }
    }
}

/// Load a set of common Restriction enzymes. Call this at program start, to load into a state field.
pub fn load_re_library() -> Vec<RestrictionEnzyme> {
    vec![
        RestrictionEnzyme::new("AatII", vec![G, A, C, G, T, C], 4),
        RestrictionEnzyme::new("Acc65I", vec![G, G, T, A, C, C], 0),
        RestrictionEnzyme::new("AscI", vec![G, G, C, G, C, G, C, C], 1),
        RestrictionEnzyme::new("AsiSI", vec![G, C, G, A, T, C, G, C], 4),
        RestrictionEnzyme::new("BamHI", vec![G, G, A, T, C, C], 0),
        RestrictionEnzyme::new("BcII", vec![T, G, A, T, C, A], 0),
        RestrictionEnzyme::new("BglII", vec![A, G, A, T, C, T], 0),
        RestrictionEnzyme::new("BmtI", vec![G, C, T, A, G, C], 4),
        RestrictionEnzyme::new("ClaI", vec![A, T, C, G, A, T], 1),
        RestrictionEnzyme::new("EcoRI", vec![G, A, A, T, T, C], 0),
        RestrictionEnzyme::new("EcoRV", vec![G, A, T, A, T, C], 2),
        RestrictionEnzyme::new("HindIII", vec![A, A, G, C, T, T], 0),
        RestrictionEnzyme::new("HpaI", vec![G, T, T, A, A, C], 2),
        RestrictionEnzyme::new("KnpI", vec![G, G, T, A, C, C], 4),
        RestrictionEnzyme::new("MscI", vec![T, G, G, C, C, A], 2),
        RestrictionEnzyme::new("NdeI", vec![C, A, T, A, T, G], 1),
        RestrictionEnzyme::new("NotI", vec![G, C, G, G, C, C, G, C], 1),
        RestrictionEnzyme::new("PsiI", vec![T, T, A, T, A, A], 2),
        RestrictionEnzyme::new("PstI", vec![C, T, G, C, A, G], 4),
        RestrictionEnzyme::new("SacI", vec![G, A, G, C, T, C], 4),
        RestrictionEnzyme::new("SalI", vec![G, T, C, G, A, C], 0),
        RestrictionEnzyme::new("SmaI", vec![C, C, C, G, G, G], 2),
        RestrictionEnzyme::new("SfbI", vec![C, C, T, G, C, A, G, G], 5),
        RestrictionEnzyme::new("SpeI", vec![A, C, T, A, G, T], 0),
        RestrictionEnzyme::new("XbaI", vec![T, C, T, A, G, A], 0),
        RestrictionEnzyme::new("XhoI", vec![C, T, C, G, A, G], 0),
        RestrictionEnzyme::new("ZraI", vec![G, A, C, G, T, C], 2),
        // RestrictionEnzyme::new("HaeIII", vec![G, G, C, C], 1), // Too many matches
    ]
}
