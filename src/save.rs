//! This module includes code for saving and loading to variou file formats.

use std::{
    fs::File,
    io,
    io::{ErrorKind, Read, Write},
};

use bincode::{config, Decode, Encode};

use crate::{primer::PrimerData, sequence::Seq, IonConcentrations, State};

pub const DEFAULT_SAVE_FILE: &str = "plasmid_tools.save";
pub const DEFAULT_FASTA_FILE: &str = "export.FASTA";

/// This is similar to `State`, but excludes the UI, and other things we don't wish to save.
#[derive(Encode, Decode)]
pub struct StateToSave {
    primer_data: Vec<PrimerData>,
    seq: Seq,
    insert_loc: usize,
    ion_concentrations: IonConcentrations,
}

impl StateToSave {
    pub fn from_state(state: &State) -> Self {
        Self {
            primer_data: state.primer_data.clone(),
            seq: state.seq.clone(),
            insert_loc: state.insert_loc,
            ion_concentrations: state.ion_concentrations.clone(),
        }
    }

    /// Used to load to state. The result is data from this struct, augmented with default values.
    pub fn to_state(self) -> State {
        State {
            primer_data: self.primer_data,
            seq: self.seq,
            insert_loc: self.insert_loc,
            ion_concentrations: self.ion_concentrations,
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
