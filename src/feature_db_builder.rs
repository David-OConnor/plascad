//! This is the entrypoint for a standalone program that parse features from GenBank and SnapGene files.
//! It assigns each sequence a label and feature type.

use std::{
    fs,
    path::{Path, PathBuf},
    str::FromStr,
};

mod file_io;
mod gui;
mod primer;
mod sequence; // Required due to naviation::Page being in file_io::save

// mod main; // todo temp
mod feature_db_load; // todo temp

use sequence::{FeatureType, Seq};

use crate::file_io::genbank;

const SOURCE_DIR: &str = "data";

struct FeatureMapItem {
    feature_name: String,
    feature_type: FeatureType,
    seq: Seq,
}

fn collect_files_in_directory(dir: &Path) -> Vec<PathBuf> {
    let mut files = Vec::new();

    if dir.is_dir() {
        for entry in fs::read_dir(dir).unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            if path.is_file() {
                files.push(path);
            }
        }
    }

    files
}

/// Combine feature maps that have identical sequences.
fn consolidate(feature_map: &[FeatureMapItem]) -> Vec<FeatureMapItem> {
    let mut result = Vec::new();

    result
}

fn main() {
    // todo: all files
    let files_genbank = collect_files_in_directory(&PathBuf::from_str(SOURCE_DIR).unwrap());
    let files_snapgene = [];

    let mut feature_maps = Vec::new();

    println!("Parsing {} files", files_genbank.len());
    for file in files_genbank {
        let path = PathBuf::from_str(file).unwrap();

        if let Ok(data) = genbank::import_genbank(&path) {
            for feature in &data.features {
                feature_maps.push({
                    FeatureMapItem {
                        name: feature.name.clone(),
                        feature_type: feature.feature_type,
                        seq: data.seq[feature.index_range.0 - 1..feature.index_range.1],
                    }
                })
            }
        }
    }

    // todo: SQlite DB?
}
