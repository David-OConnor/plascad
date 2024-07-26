//! This module contains info related to Restriction Enzyme sites.
//!
//! [Wikipedia: List of RE sites](https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites:_A)
//! [NEB guide](https://www.neb.com/en-us/tools-and-resources/selection-charts/frequencies-of-restriction-sites)

use bincode::{Decode, Encode};

use crate::{
    primer::PrimerDirection,
    Nucleotide::{A, C, G, T},
    Seq,
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
}

#[derive(Encode, Decode)]
pub struct ReMatch {
    pub lib_index: usize,
    pub seq_index: usize,
    /// Direction helps align the cut site.
    pub direction: PrimerDirection,
}

#[derive(Encode, Decode)] // todo: IDeally, remove this.
pub struct RestrictionEnzyme {
    pub name: String,
    /// From the 5' end.
    pub seq: Seq, // todo: General?
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

pub fn get_common_res() -> Vec<RestrictionEnzyme> {
    vec![
        RestrictionEnzyme::new("AatII", vec![G, A, C, G, T, C], 4),
        RestrictionEnzyme::new("BamHI", vec![G, G, A, T, C, C], 0),
        RestrictionEnzyme::new("BcII", vec![T, G, A, T, C, A], 0),
        RestrictionEnzyme::new("EcoRI", vec![G, A, C, G, T, C], 0),
        RestrictionEnzyme::new("EcoRV", vec![G, A, T, A, T, C], 2),
        RestrictionEnzyme::new("HindIII", vec![A, A, G, C, T, T], 0),
        RestrictionEnzyme::new("HpaI", vec![G, T, T, A, A, C], 2),
        RestrictionEnzyme::new("KnpI", vec![G, G, T, A, C, C], 4),
        RestrictionEnzyme::new("MscI", vec![T, G, G, C, C, A], 2),
        RestrictionEnzyme::new("NdeI", vec![C, A, T, A, T, G], 1),
        RestrictionEnzyme::new("NotI", vec![G, C, G, G, C, C, G, C], 1),
        RestrictionEnzyme::new("PstI", vec![C, T, G, C, A, G], 4),
        RestrictionEnzyme::new("SaII", vec![G, T, C, G, A, C], 0),
        RestrictionEnzyme::new("SmaI", vec![C, C, C, G, G, G], 2),
        RestrictionEnzyme::new("SpeI", vec![A, C, T, A, G, T], 0),
        RestrictionEnzyme::new("XhoI", vec![C, T, C, G, A, G], 0),
        RestrictionEnzyme::new("HaeIII", vec![G, G, C, C], 1),
    ]
}
