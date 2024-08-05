//! This module contains code for saving and loading in several file formats.

use crate::{
    primer::Primer,
    sequence::{Feature, Seq, SeqTopology},
    Reference,
};

pub mod genbank;
pub mod save;
pub mod snapgene;

/// The most important data to store, used by our format, GenBank, and SnapGene.
/// todo: Integrate this directly into our State struct, A/R.
#[derive(Default)]
pub struct GenericData {
    pub seq: Seq,
    pub plasmid_name: String,
    pub topology: SeqTopology,
    pub features: Vec<Feature>,
    pub primers: Vec<Primer>,
    pub comments: Vec<String>,
    pub references: Vec<Reference>,
}
