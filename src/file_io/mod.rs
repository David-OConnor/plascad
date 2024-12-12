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
        .add_file_filter(
            "PlasCAD files",
            Arc::new(|p| p.extension().unwrap_or_default().to_ascii_lowercase() == "pcad"),
        )
        .add_file_filter(
            "FASTA files",
            Arc::new(|p| p.extension().unwrap_or_default().to_ascii_lowercase() == "fasta"),
        )
        .add_file_filter(
            "GenBank files",
            Arc::new(|p| {
                let ext = p.extension().unwrap_or_default().to_ascii_lowercase();
                ext == "gb" || ext == "gbk"
            }),
        )
        .add_file_filter(
            "SnapGene DNA files",
            Arc::new(|p| p.extension().unwrap_or_default().to_ascii_lowercase() == "dna"),
        )
        .add_file_filter(
            // Note: We experience glitches if this name is too long. (Window extends horizontally)
            "PCAD/FASTA/GB/SG/AB1",
            Arc::new(|p| {
                let ext = p.extension().unwrap_or_default().to_ascii_lowercase();
                ext == "pcad"
                    || ext == "fasta"
                    || ext == "gb"
                    || ext == "gbk"
                    || ext == "dna"
                    || ext == "ab1"
            }),
        );

        let save = FileDialog::new()
            // .add_quick_access("Project", |s| {
            //     s.add_path("☆  Examples", "examples");
            // })
            .add_file_filter(
                "PlasCAD files",
                Arc::new(|p| p.extension().unwrap_or_default().to_ascii_lowercase() == "pcad"),
            )
            .default_file_filter("PlasCAD files")
            .default_file_name(QUICKSAVE_FILE)
            .id("0");

        let import = FileDialog::with_config(cfg_import.clone())
            .default_file_filter("PCAD/FASTA/GB/SG")
            .id("1");

        let export_fasta = FileDialog::new()
            .add_file_filter(
                "FASTA files",
                Arc::new(|p| p.extension().unwrap_or_default().to_ascii_lowercase() == "fasta"),
            )
            .default_file_filter("FASTA files")
            .default_file_name(DEFAULT_FASTA_FILE)
            .id("3");

        let export_genbank = FileDialog::new()
            .add_file_filter(
                "GenBank files",
                Arc::new(|p| {
                    let ext = p.extension().unwrap_or_default().to_ascii_lowercase();
                    ext == "gb" || ext == "gbk"
                }),
            )
            .default_file_filter("GenBank files")
            .default_file_name(DEFAULT_GENBANK_FILE)
            .id("4");

        let export_dna = FileDialog::new()
            .add_file_filter(
                "SnapGene DNA files",
                Arc::new(|p| p.extension().unwrap_or_default().to_ascii_lowercase() == "dna"),
            )
            .default_file_filter("SnapGene DNA files")
            .default_file_name(DEFAULT_DNA_FILE)
            .id("5");

        let cloning_import = FileDialog::with_config(cfg_import)
            .default_file_filter("PCAD/FASTA/GB/SG")
            .id("6");

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
