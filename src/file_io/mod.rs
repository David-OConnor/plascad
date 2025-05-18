//! This module contains code for saving and loading in several file formats.

use std::{path::Path, sync::Arc};

use egui_file_dialog::{FileDialog, FileDialogConfig};
use na_seq::{Seq, SeqTopology};

use crate::{
    file_io::save::{DEFAULT_DNA_FILE, DEFAULT_FASTA_FILE, DEFAULT_GENBANK_FILE, QUICKSAVE_FILE},
    misc_types::{Feature, Metadata},
    primer::Primer,
};

pub mod ab1;
mod ab1_tags;
pub mod genbank;
mod pcad;
pub mod save;
pub mod snapgene;

/// The most important data to store, used by our format, GenBank, and SnapGene.
/// We use this in our main State struct to keep track of this data.
#[derive(Default, Clone)]
pub struct GenericData {
    pub seq: Seq,
    pub topology: SeqTopology,
    pub features: Vec<Feature>,
    pub primers: Vec<Primer>,
    pub metadata: Metadata,
}

pub struct FileDialogs {
    pub save: FileDialog,
    pub load: FileDialog,
    pub export_fasta: FileDialog,
    pub export_genbank: FileDialog,
    pub export_dna: FileDialog,
    pub cloning_load: FileDialog,
}

impl Default for FileDialogs {
    fn default() -> Self {
        // We can't clone `FileDialog`; use this to reduce repetition instead.
        let cfg_import = FileDialogConfig {
            // todo: Explore other optiosn A/R
            ..Default::default()
        }
        .add_file_filter_extensions("PlasCAD", vec!["pcad"])
        .add_file_filter_extensions("FASTA", vec!["fasta"])
        .add_file_filter_extensions("GenBank", vec!["gb, gbk"])
        .add_file_filter_extensions("SnapGene DNA", vec!["dna"])
        .add_file_filter_extensions(
            "PCAD/FASTA/GB/DNA/AB1",
            vec!["pcad", "fasta, gb, gbk, dna, ab1"],
        );

        let save = FileDialog::new()
            // .add_quick_access("Project", |s| {
            //     s.add_path("☆  Examples", "examples");
            // })
            .add_save_extension("PlasCAD", "pcad")
            .default_save_extension("PlasCAD")
            .default_file_name(QUICKSAVE_FILE);
        // .id("0");

        let import = FileDialog::with_config(cfg_import.clone())
            .default_file_filter("PCAD/FASTA/GB/DNA/AB1");

        let export_fasta = FileDialog::new()
            .add_save_extension("FASTA", "fasta")
            .default_save_extension("FASTA")
            .default_file_name(DEFAULT_FASTA_FILE);

        let export_genbank = FileDialog::new()
            .add_save_extension("GenBank", "gb")
            .default_save_extension("GenBank")
            .default_file_name(DEFAULT_GENBANK_FILE);

        let export_dna = FileDialog::new()
            .add_save_extension("SnapGene DNA", "dna")
            .default_save_extension("SnapGene DNA")
            .default_file_name(DEFAULT_DNA_FILE);

        let cloning_import =
            FileDialog::with_config(cfg_import).default_file_filter("PCAD/FASTA/GB/DNA/AB1");

        Self {
            save,
            // load: load_,
            load: import,
            export_fasta,
            export_genbank,
            export_dna,
            cloning_load: cloning_import,
            // selected: None,
        }
    }
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
