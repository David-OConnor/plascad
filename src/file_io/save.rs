//! This module includes code for reading and writing our PlasCAD binary fileformat.
//!
//! Todo: Break into packets. This may require further ejecting from bincode, at least at the top level[s].

use std::{
    fs::{File, OpenOptions},
    io,
    io::{ErrorKind, Read, Write},
    path::Path,
};

use bincode::{
    config,
    error::{DecodeError, EncodeError},
    BorrowDecode, Decode, Encode,
};
use bio::io::fasta;

use crate::{
    file_io::GenericData,
    gui::navigation::{Page, PageSeq, PageSeqTop},
    primer::Primer,
    sequence::{Feature, Nucleotide, ReadingFrame, Seq, SeqTopology},
    IonConcentrations, Metadata, PcrUi, SeqVisibility, State, StateUi,
};

pub const DEFAULT_SAVE_FILE: &str = "plasmid.pcad";
pub const DEFAULT_PREFS_FILE: &str = "pcad_prefs.pp";

pub const DEFAULT_FASTA_FILE: &str = "export.fasta";
pub const DEFAULT_GENBANK_FILE: &str = "export.gbk";
pub const DEFAULT_DNA_FILE: &str = "export.dna";

/// This is similar to `State`, but excludes the UI, and other things we don't wish to save.
pub struct StateToSave {
    generic: GenericData,
    insert_loc: usize,
    ion_concentrations: IonConcentrations,
    reading_frame: ReadingFrame,
}

impl Encode for StateToSave {
    fn encode<E: bincode::enc::Encoder>(&self, encoder: &mut E) -> Result<(), EncodeError> {
        // todo: We currently use default encoding for all this, but may change later.
        self.generic.encode(encoder)?;
        self.insert_loc.encode(encoder)?;
        self.ion_concentrations.encode(encoder)?;
        self.reading_frame.encode(encoder)?;

        Ok(())
    }
}

impl Decode for StateToSave {
    fn decode<D: bincode::de::Decoder>(decoder: &mut D) -> Result<Self, DecodeError> {
        // todo: We currently use default encoding for all this, but may change later.
        let generic = GenericData::decode(decoder)?;
        let insert_loc = usize::decode(decoder)?;
        let ion_concentrations = IonConcentrations::decode(decoder)?;
        let reading_frame = ReadingFrame::decode(decoder)?;

        Ok(Self {
            generic,
            insert_loc,
            ion_concentrations,
            reading_frame,
        })
    }
}

impl Encode for GenericData {
    fn encode<E: bincode::enc::Encoder>(&self, encoder: &mut E) -> Result<(), EncodeError> {
        // Serialize seq using our custom serializer
        let seq_data = serialize_seq_bin(&self.seq);
        seq_data.encode(encoder)?;

        // Serialize other fields using default serialization
        self.topology.encode(encoder)?;
        self.features.encode(encoder)?;
        self.primers.encode(encoder)?;
        self.metadata.encode(encoder)?;

        Ok(())
    }
}

impl Decode for GenericData {
    fn decode<D: bincode::de::Decoder>(decoder: &mut D) -> Result<Self, DecodeError> {
        // Deserialize seq using our custom deserializer
        let seq_data = Vec::decode(decoder)?;
        // let seq = deser_seq_bin(&seq_data).map_err(|_| file_io: :Error::new(ErrorKind::InvalidData, "Error deserializing seq"));
        let seq = deser_seq_bin(&seq_data).unwrap_or_default(); // todo: Better error handling/propogation.

        // Deserialize other fields using default deserialization
        let topology = SeqTopology::decode(decoder)?;
        let features = Vec::<Feature>::decode(decoder)?;
        let primers = Vec::<Primer>::decode(decoder)?;
        let metadata = Metadata::decode(decoder)?;

        Ok(Self {
            seq,
            topology,
            features,
            primers,
            metadata,
        })
    }
}

impl StateToSave {
    pub fn from_state(state: &State) -> Self {
        Self {
            generic: state.generic.clone(),
            insert_loc: state.insert_loc,
            ion_concentrations: state.ion_concentrations.clone(),
            reading_frame: state.reading_frame.clone(),
        }
    }

    /// Used to load to state. The result is data from this struct, augmented with default values.
    pub fn to_state(self) -> State {
        State {
            generic: self.generic,
            insert_loc: self.insert_loc,
            ion_concentrations: self.ion_concentrations,
            reading_frame: self.reading_frame,
            ..Default::default()
        }
    }
}

#[derive(Encode, Decode)]
pub struct StateUiToSave {
    page: Page,
    page_seq: PageSeq,
    page_seq_top: PageSeqTop,
    pcr: PcrUi,
    primer_selected: Option<usize>,
    feature_selected: Option<usize>,
    seq_visibility: SeqVisibility,
    hide_map_feature_editor: bool,
}

impl StateUiToSave {
    pub fn from_state(state: &StateUi) -> Self {
        Self {
            page: state.page.clone(),
            page_seq: state.page_seq.clone(),
            page_seq_top: state.page_seq_top.clone(),
            pcr: state.pcr.clone(),
            primer_selected: state.primer_selected.clone(),
            feature_selected: state.feature_selected.clone(),
            seq_visibility: state.seq_visibility.clone(),
            hide_map_feature_editor: state.hide_map_feature_editor.clone(),
        }
    }

    /// Used to load to state. The result is data from this struct, augmented with default values.
    pub fn to_state(self) -> StateUi {
        StateUi {
            page: self.page,
            page_seq: self.page_seq,
            page_seq_top: self.page_seq_top,
            pcr: self.pcr,
            primer_selected: self.primer_selected,
            feature_selected: self.feature_selected,
            seq_visibility: self.seq_visibility,
            hide_map_feature_editor: self.hide_map_feature_editor,
            ..Default::default()
        }
    }
}

pub fn save<T: Encode>(path: &Path, data: &T) -> io::Result<()> {
    let config = config::standard();

    let encoded: Vec<u8> = bincode::encode_to_vec(data, config).unwrap();

    let mut file = File::create(path)?;
    file.write_all(&encoded)?;
    Ok(())
}

pub fn load<T: Decode>(path: &Path) -> io::Result<T> {
    let config = config::standard();

    let mut file = File::open(path)?;
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
    let file = File::create(path)?;

    let seq_u8: Vec<u8> = seq.iter().map(|nt| nt.to_u8_letter()).collect();
    let mut writer = fasta::Writer::new(file);

    writer.write(name, Some("A DNA export from PlasCAD"), seq_u8.as_slice())?;

    Ok(())
}

/// Import from a FASTA file. (Seq, plasmid name (id), description)
pub fn import_fasta(path: &Path) -> io::Result<(Seq, String, String)> {
    let file = File::open(path)?;

    let mut records = fasta::Reader::new(file).records();

    let mut result = Vec::new();

    // todo: Do we want id, or description?
    let mut id = String::new();
    let mut description = String::new();

    while let Some(Ok(record)) = records.next() {
        for r in record.seq() {
            result.push(Nucleotide::from_u8_letter(*r)?);
            record.id().clone_into(&mut id); // Note that this overrides previous records, if applicable.
            record
                .desc()
                .unwrap_or_default()
                .clone_into(&mut description)
        }
    }

    Ok((result, id, description))
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
