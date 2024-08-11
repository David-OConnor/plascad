//! This module contains code for saving and loading in several file formats.

use std::path::Path;

use crate::{
    primer::Primer,
    sequence::{Feature, Metadata, Seq, SeqTopology},
};

pub mod genbank;
pub mod save;
pub mod snapgene;

/// The most important data to store, used by our format, GenBank, and SnapGene.
/// todo: Integrate this directly into our State struct, A/R.
#[derive(Default, Clone)]
pub struct GenericData {
    pub seq: Seq,
    pub topology: SeqTopology,
    pub features: Vec<Feature>,
    pub primers: Vec<Primer>,
    pub metadata: Metadata,
}

/// There doesn't seem to be a clear name in GenBank or Snapgene formats; use the filename.
/// Note: This includes error checking, but this should always pass under normal circumstances.
fn get_filename(path: &Path) -> String {
    if let Some(file_name) = path.file_stem() {
        file_name
            .to_str()
            .map(|s| s.to_string())
            .unwrap_or_default()
    } else {
        String::new()
    }
}
