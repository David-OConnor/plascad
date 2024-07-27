//! This module includes code for saving and loading to variou file formats.

use std::{
    fs::{File, OpenOptions},
    io,
    io::{ErrorKind, Read, Write},
    path::Path,
};

use bincode::{config, Decode, Encode};
use bio::io::fasta;

use crate::{
    primer::PrimerData,
    sequence::{Feature, Nucleotide, Seq},
    IonConcentrations, State,
};

pub const DEFAULT_SAVE_FILE: &str = "plasmid_tools.save";
pub const DEFAULT_FASTA_FILE: &str = "export.fasta";

/// This is similar to `State`, but excludes the UI, and other things we don't wish to save.
#[derive(Encode, Decode)]
pub struct StateToSave {
    primer_data: Vec<PrimerData>,
    seq: Seq,
    insert_loc: usize,
    ion_concentrations: IonConcentrations,
    features: Vec<Feature>,
}

impl StateToSave {
    pub fn from_state(state: &State) -> Self {
        Self {
            primer_data: state.primer_data.clone(),
            seq: state.seq.clone(),
            insert_loc: state.insert_loc,
            ion_concentrations: state.ion_concentrations.clone(),
            features: state.features.clone(),
        }
    }

    /// Used to load to state. The result is data from this struct, augmented with default values.
    pub fn to_state(self) -> State {
        State {
            primer_data: self.primer_data,
            seq: self.seq,
            insert_loc: self.insert_loc,
            ion_concentrations: self.ion_concentrations,
            features: self.features,
            ..Default::default()
        }
    }
}

pub fn save<T: Encode>(filename: &str, data: &T) -> io::Result<()> {
    let config = config::standard();

    let encoded: Vec<u8> = bincode::encode_to_vec(data, config).unwrap();
    let mut file = File::create(filename)?;
    file.write_all(&encoded)?;
    Ok(())
}

pub fn load<T: Decode>(filename: &str) -> io::Result<T> {
    let config = config::standard();

    let mut file = File::open(filename)?;
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;
    let (decoded, len) = match bincode::decode_from_slice(&buffer, config) {
        Ok(v) => v,
        Err(_) => {
            eprintln!("Error loading from file. Did the format change?");
            return Err(io::Error::new(ErrorKind::Other, "error loading"));
        }
    };
    Ok(decoded)
}

/// Export a sequence in FASTA format.
pub fn export_fasta(seq: &[Nucleotide], path: &Path) -> io::Result<()> {
    let file = OpenOptions::new()
        .append(true)
        .create(true) // Create the file if it doesn't exist
        .open(path)?;

    let seq_u8: Vec<u8> = seq.iter().map(|nt| *nt as u8).collect();
    let mut writer = fasta::Writer::new(file);

    writer.write(
        "Plasmid tools export",
        Some("A DNA export from Plasmid tools"),
        seq_u8.as_slice(),
    )?;

    Ok(())
}

/// Import from a FASTA file
pub fn import_fasta(path: &Path) -> io::Result<Seq> {
    let file = File::open(path)?;

    let mut records = fasta::Reader::new(file).records();

    let mut result = Vec::new();

    while let Some(Ok(record)) = records.next() {
        for r in record.seq() {
            result.push(Nucleotide::from_u8(*r)?);
        }
    }

    Ok(result)
}
