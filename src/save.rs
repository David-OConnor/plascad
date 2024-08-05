//! This module includes code for saving and loading to variou file formats.

use std::{
    fs::{File, OpenOptions},
    io,
    io::{ErrorKind, Read, Write},
    path::Path,
};

use bincode::{
    config,
    error::{DecodeError, EncodeError},
    Decode, Encode,
};
use bio::io::fasta;

use crate::{primer::PrimerData, sequence::{Feature, Nucleotide, ReadingFrame, Seq, SeqTopology}, IonConcentrations, State, Reference};

pub const DEFAULT_SAVE_FILE: &str = "plasmid_tools.pcad";
pub const DEFAULT_FASTA_FILE: &str = "export.fasta";
pub const DEFAULT_DNA_FILE: &str = "export.dna";

/// This is similar to `State`, but excludes the UI, and other things we don't wish to save.
// #[derive(Encode, Decode)]
pub struct StateToSave {
    seq: Seq,
    features: Vec<Feature>,
    primer_data: Vec<PrimerData>,
    insert_loc: usize,
    ion_concentrations: IonConcentrations,
    plasmid_name: String,
    topology: SeqTopology,
    reading_frame: ReadingFrame,
    comments: Vec<String>,
    references: Vec<Reference>,
}

impl Encode for StateToSave {
    fn encode<E: bincode::enc::Encoder>(&self, encoder: &mut E) -> Result<(), EncodeError> {
        // Serialize seq using our custom serializer
        let seq_data = serialize_seq_bin(&self.seq);
        seq_data.encode(encoder)?;

        // Serialize other fields using default serialization
        self.features.encode(encoder)?;
        self.primer_data.encode(encoder)?;
        self.insert_loc.encode(encoder)?;
        self.ion_concentrations.encode(encoder)?;
        self.plasmid_name.encode(encoder)?;
        self.topology.encode(encoder)?;
        self.reading_frame.encode(encoder)?;
        self.comments.encode(encoder)?;
        self.references.encode(encoder)?;

        Ok(())
    }
}

impl Decode for StateToSave {
    fn decode<D: bincode::de::Decoder>(decoder: &mut D) -> Result<Self, DecodeError> {
        // Deserialize seq using our custom deserializer
        let seq_data = Vec::decode(decoder)?;
        // let seq = deser_seq_bin(&seq_data).map_err(|_| io: :Error::new(ErrorKind::InvalidData, "Error deserializing seq"));
        let seq = deser_seq_bin(&seq_data).unwrap_or_default(); // todo: Better error handling/propogation.

        // Deserialize other fields using default deserialization
        let features = Vec::<Feature>::decode(decoder)?;
        let primer_data = Vec::<PrimerData>::decode(decoder)?;
        let insert_loc = usize::decode(decoder)?;
        let ion_concentrations = IonConcentrations::decode(decoder)?;
        let plasmid_name = String::decode(decoder)?;
        let topology = SeqTopology::decode(decoder)?;
        let reading_frame = ReadingFrame::decode(decoder)?;
        let comments = Vec::<String>::decode(decoder)?;
        let references = Vec::<Reference>::decode(decoder)?;

        Ok(StateToSave {
            seq,
            features,
            primer_data,
            insert_loc,
            ion_concentrations,
            plasmid_name,
            topology,
            reading_frame,
            comments,
            references,
        })
    }
}

impl StateToSave {
    pub fn from_state(state: &State) -> Self {
        Self {
            seq: state.seq.clone(),
            features: state.features.clone(),
            primer_data: state.primer_data.clone(),
            insert_loc: state.insert_loc,
            ion_concentrations: state.ion_concentrations.clone(),
            plasmid_name: state.plasmid_name.clone(),
            topology: state.topology.clone(),
            reading_frame: state.reading_frame.clone(),
            comments: state.comments.clone(),
            references: state.references.clone(),
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
            plasmid_name: self.plasmid_name,
            topology: self.topology,
            reading_frame: self.reading_frame,
            comments: self.comments,
            references: self.references,
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
    let (decoded, _len) = match bincode::decode_from_slice(&buffer, config) {
        Ok(v) => v,
        Err(_) => {
            eprintln!("Error loading from file. Did the format change?");
            return Err(io::Error::new(ErrorKind::Other, "error loading"));
        }
    };
    Ok(decoded)
}

/// Export a sequence in FASTA format.
pub fn export_fasta(seq: &[Nucleotide], name: &str, path: &Path) -> io::Result<()> {
    let file = OpenOptions::new()
        .write(true)
        .create(true) // Create the file if it doesn't exist
        .open(path)?;

    let seq_u8: Vec<u8> = seq.iter().map(|nt| nt.to_u8_letter()).collect();
    let mut writer = fasta::Writer::new(file);

    writer.write(name, Some("A DNA export from PlasCAD"), seq_u8.as_slice())?;

    Ok(())
}

/// Import from a FASTA file
pub fn import_fasta(path: &Path) -> io::Result<(Seq, String)> {
    let file = File::open(path)?;

    let mut records = fasta::Reader::new(file).records();

    let mut result = Vec::new();

    // todo: Do we want id, or description?
    let mut id = String::new();

    while let Some(Ok(record)) = records.next() {
        for r in record.seq() {
            result.push(Nucleotide::from_u8_letter(*r)?);
            id = record.id().to_owned(); // Note that this overrides previous records, if applicable.
        }
    }

    Ok((result, id))
}

/// A compact binary serialization of our sequence. Useful for file storage.
/// The first byte is sequence length; we need this, since one of our nucleotides necessarily serializes
/// to 0b00.
/// todo: Is this MSB or LSB?
pub fn serialize_seq_bin(seq: &[Nucleotide]) -> Vec<u8> {
    let mut result = Vec::new();
    result.extend(&(seq.len() as u32).to_be_bytes());

    for i in 0..seq.len() / 4 + 1 {
        let mut val = 0;
        for j in 0..4 {
            let ind = i * 4 + j;
            if ind + 1 > seq.len() {
                break;
            }
            let nt = seq[ind];
            val |= (nt as u8) << (j * 2);
        }
        result.push(val);
    }
    result
}

/// A compact binary deserialization of our sequence.
/// todo: Is this MSB or LSB?
pub fn deser_seq_bin(data: &[u8]) -> io::Result<Seq> {
    let mut result = Vec::new();

    if data.len() < 4 {
        return Err(io::Error::new(
            ErrorKind::InvalidData,
            "Bin nucleotide sequence is too short.",
        ));
    }

    let seq_len = u32::from_be_bytes(data[0..4].try_into().unwrap()) as usize;

    for byte in &data[4..] {
        for i in 0..4 {
            // This trimming removes extra 00-serialized nucleotides.
            if result.len() >= seq_len {
                break;
            }

            let bits = (byte >> (2 * i)) & 0b11;
            result.push(Nucleotide::try_from(bits).map_err(|_| {
                io::Error::new(
                    ErrorKind::InvalidData,
                    format!("Invalid NT serialization: {}, {}", byte, bits),
                )
            })?);
        }
    }

    Ok(result)
}
