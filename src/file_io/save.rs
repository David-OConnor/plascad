//! This module includes code for reading and writing our PlasCAD binary fileformat.
//!
//! Todo: Break into packets. This may require further ejecting from bincode, at least at the top level[s].

use std::{
    env,
    fs::File,
    io,
    io::{ErrorKind, Read, Write},
    path::{Path, PathBuf},
};

use bincode::{
    config,
    error::{DecodeError, EncodeError},
    Decode, Encode,
};
use bio::io::fasta;
use eframe::egui::Ui;
use na_seq::{deser_seq_bin, seq_to_u8_lower, serialize_seq_bin, Nucleotide, Seq, SeqTopology};

use crate::{
    ab1::SeqRecordAb1,
    feature_db_load::find_features,
    file_io::{
        ab1::import_ab1,
        genbank::{export_genbank, import_genbank},
        snapgene::{export_snapgene, import_snapgene},
        GenericData,
    },
    gui::{
        navigation::{Page, PageSeq, PageSeqTop, Tab},
        set_window_title,
    },
    misc_types::{Feature, Metadata},
    pcr::PcrUi,
    portions::PortionsState,
    primer::{IonConcentrations, Primer},
    state::State,
    Selection, SeqVisibility, StateUi,
};

pub const QUICKSAVE_FILE: &str = "quicksave.pcad";
pub const DEFAULT_PREFS_FILE: &str = "pcad_prefs.pp";

pub const DEFAULT_FASTA_FILE: &str = "export.fasta";
pub const DEFAULT_GENBANK_FILE: &str = "export.gbk";
pub const DEFAULT_DNA_FILE: &str = "export.dna";

/// Sequence-related data to save in our own file format, GBK, or Snapgene.
#[derive(Default)]
pub struct StateToSave {
    pub generic: GenericData,
    // pub path_loaded: Option<Tab>,
    // pub ion_concentrations: IonConcentrations,
    pub portions: PortionsState,
    // pub ab1_data: Vec<SeqRecordAb1>,
    // todo: Sort this out; how to indicate we loaded AB1 vs normal.
    pub ab1_data: SeqRecordAb1,
    pub path_loaded: Option<PathBuf>,
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
    pub fn from_state(state: &State, active: usize) -> Self {
        Self {
            generic: state.generic[active].clone(),
            // insert_loc: state.cloning_insert_loc, // todo: Not fully handled.
            // ion_concentrations: state.ion_concentrations[active].clone(),
            // path_loaded: state.path_loaded[active].clone(),
            portions: state.portions[state.active].clone(),
            // ab1_data: state.ab1_data.clone(),
            ab1_data: state.ab1_data[state.active].clone(),
            path_loaded: None, // todo: Is this correct?
        }
    }

    /// Saves in PCAD format. todo: Move to the PCAD file A/R.
    pub fn save_to_file(&self, path: &Path) -> io::Result<()> {
        let encoded = self.to_bytes();

        let mut file = File::create(path)?;
        file.write_all(&encoded)?;
        Ok(())
    }

    /// Loads in PCAD format.
    pub fn load_from_file(path: &Path) -> io::Result<Self> {
        let mut file = File::open(path)?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;

        Self::from_bytes(&buffer)
    }
}

#[derive(Encode, Decode)]
/// Represents state to save automatically; not related to sequence data.
pub struct PrefsToSave {
    page: Page,
    page_seq: PageSeq,
    page_seq_top: PageSeqTop,
    pcr: PcrUi,
    selected_item: Selection,
    seq_visibility: SeqVisibility,
    hide_map_feature_editor: bool,
    tabs_open: Vec<Tab>,
    ion_concentrations: IonConcentrations,
}

impl PrefsToSave {
    pub fn from_state(
        state: &StateUi,
        tabs_open_: &[Tab],
        ion_concentrations: &IonConcentrations,
    ) -> Self {
        // Remove the empty paths; we can't load them.
        let mut tabs_open = Vec::new();
        for t in tabs_open_ {
            if t.path.is_some() {
                tabs_open.push(t.clone());
            }
        }

        Self {
            page: state.page,
            page_seq: state.page_seq,
            page_seq_top: state.page_seq_top,
            pcr: state.pcr.clone(),
            selected_item: state.selected_item,
            seq_visibility: state.seq_visibility.clone(),
            hide_map_feature_editor: state.hide_map_feature_editor,
            tabs_open,
            ion_concentrations: ion_concentrations.clone(),
        }
    }

    /// Used to load to state. The result is data from this struct, augmented with default values.
    pub fn to_state(&self) -> (StateUi, Vec<Tab>, IonConcentrations) {
        (
            StateUi {
                page: self.page,
                page_seq: self.page_seq,
                page_seq_top: self.page_seq_top,
                pcr: self.pcr.clone(),
                selected_item: self.selected_item,
                seq_visibility: self.seq_visibility.clone(),
                hide_map_feature_editor: self.hide_map_feature_editor,
                // last_file_opened: self.last_file_opened.clone(),
                ..Default::default()
            },
            self.tabs_open.clone(),
            self.ion_concentrations.clone(),
        )
    }
}

/// Save to file, using Bincode. We currently use this for preference files.
pub fn save<T: Encode>(path: &Path, data: &T) -> io::Result<()> {
    let config = config::standard();

    let encoded: Vec<u8> = bincode::encode_to_vec(data, config).unwrap();

    let mut file = File::create(path)?;
    file.write_all(&encoded)?;
    Ok(())
}

/// Load from file, using Bincode. We currently use this for preference files.
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

    let mut writer = fasta::Writer::new(file);

    writer.write(
        name,
        Some("A DNA export from PlasCAD"),
        seq_to_u8_lower(seq).as_slice(),
    )?;

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

/// Save a new file, eg a cloning or PCR product.
pub fn save_new_product(name: &str, state: &mut State, ui: &mut Ui) {
    name.clone_into(&mut state.generic[state.active].metadata.plasmid_name);

    let active = state.generic.len();
    state.active = active;

    // todo: Option for GenBank and SnapGene formats here?
    let mut save_path = env::current_dir().unwrap();
    let filename = {
        let name = state.generic[state.active]
            .metadata
            .plasmid_name
            .to_lowercase()
            .replace(' ', "_");
        format!("{name}.pcad")
    };
    save_path.push(Path::new(&filename));

    state.ui.file_dialogs.save.config_mut().default_file_name = filename.to_string();
    state.ui.file_dialogs.save.save_file();

    if let Some(path) = state.ui.file_dialogs.save.take_selected() {
        match StateToSave::from_state(state, state.active).save_to_file(&path) {
            Ok(_) => {
                // state.file_active = Some(Tab {
                //     path: path.to_owned(),
                //     ab1: false,
                // });

                set_window_title(&state.tabs_open[state.active], ui);

                // if let Some(path) = state.tabs_open[state.active].path {
                //
                // }
            }
            Err(e) => eprintln!(
                "Error saving cloning product in the PlasCAD format: {:?}",
                e
            ),
        };
    }
}

/// Load state from a file of various formats.
pub fn load_import(path: &Path) -> Option<StateToSave> {
    let mut result = StateToSave::default();

    if let Some(extension) = path.extension().and_then(|ext| ext.to_str()) {
        match extension.to_lowercase().as_ref() {
            "pcad" => {
                let state_loaded = StateToSave::load_from_file(path);
                match state_loaded {
                    Ok(s) => {
                        result.generic = s.generic;
                        result.path_loaded = Some(path.to_owned());
                        // result.ion_concentrations = s.ion_concentrations;
                        // result.path_loaded = Some(Tab {
                        //     path: path.to_owned(),
                        //     ab1: false,
                        // });
                        result.portions = s.portions;

                        return Some(result);
                    }
                    Err(e) => {
                        eprintln!("Error loading a PCAD file: {e}");
                    }
                };
            }
            // Does this work for FASTQ too?
            "fasta" | "fa" => {
                if let Ok((seq, id, description)) = import_fasta(path) {
                    result.generic.seq = seq;
                    result.generic.metadata.plasmid_name = id;
                    result.generic.metadata.comments = vec![description];
                    // FASTA is seq-only data, so don't attempt to save over it.

                    // Automatically annotate FASTA files.
                    result.generic.features = find_features(&result.generic.seq);

                    return Some(result);
                }
            }
            "dna" => {
                if let Ok(data) = import_snapgene(path) {
                    result.generic = data;
                    // We do not mark the path as opened if using SnapGene, since we currently can not
                    // fully understand the format, nor make a native file SnapGene can open.

                    return Some(result);
                }
            }
            "gb" | "gbk" => {
                if let Ok(data) = import_genbank(path) {
                    result.generic = data;
                    // result.path_loaded = Some(Tab {
                    //     path: path.to_owned(),
                    //     ab1: false,
                    // });
                    result.path_loaded = Some(path.to_owned());
                    return Some(result);
                }
            }
            "ab1" => {
                if let Ok(data) = import_ab1(path) {
                    if data.len() >= 1 {
                        result.ab1_data = data[0].clone(); // todo: Note that this assumes len 1 of results.
                    }
                    // result.path_loaded = Some(Tab {
                    //     path: path.to_owned(),
                    //     ab1: true,
                    // });
                    result.path_loaded = Some(path.to_owned());
                    return Some(result);
                }
            }

            _ => {
                eprintln!(
                    "The file to import must be in PlasCAD, FASTA, GenBank, SnapGene, or AB1 format."
                )
            }
        }
    }
    None
}

/// Save the current file ("save" vice "save as") if there is one; if not, quicksave to an anonymous file.
pub fn save_current_file(state: &State) {
    match &state.tabs_open[state.active].path {
        Some(path) => {
            if let Some(extension) = path.extension().and_then(|ext| ext.to_str()) {
                match extension.to_lowercase().as_ref() {
                    "pcad" => {
                        // if let Err(e) = save(path, &StateToSave::from_state(state, state.active)) {
                        if let Err(e) =
                            StateToSave::from_state(state, state.active).save_to_file(&path)
                        {
                            eprintln!("Error saving in PlasCAD format: {}", e);
                        };
                    }
                    // Does this work for FASTQ too?
                    "fasta" => {
                        if let Err(e) = export_fasta(
                            state.get_seq(),
                            &state.generic[state.active].metadata.plasmid_name,
                            &path,
                        ) {
                            eprintln!("Error exporting to FASTA: {:?}", e);
                        };
                    }
                    "dna" => {
                        if let Err(e) = export_snapgene(&state.generic[state.active], path) {
                            eprintln!("Error exporting to SnapGene: {:?}", e);
                        };
                    }
                    "gb" | "gbk" => {
                        let mut primer_matches = Vec::new();
                        for primer in &state.generic[state.active].primers {
                            for prim_match in &primer.volatile.matches {
                                primer_matches.push((prim_match.clone(), primer.name.clone()));
                            }
                        }

                        if let Err(e) =
                            export_genbank(&state.generic[state.active], &primer_matches, path)
                        {
                            eprintln!("Error exporting to GenBank: {:?}", e);
                        };
                    }
                    _ => {
                        eprintln!("Unexpected file format loading.")
                    }
                }
            }
        }
        None => {
            // Quicksave.
            if let Err(e) = StateToSave::from_state(state, state.active)
                .save_to_file(&PathBuf::from(QUICKSAVE_FILE))
            {
                eprintln!("Error quicksaving: {e}");
            }
        }
    }
}
