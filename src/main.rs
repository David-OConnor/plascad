// Disables the terminal window. Use this for releases, but disable when debugging.
// #![windows_subsystem = "windows"]

use std::{io, path::PathBuf, sync::Arc};

use bincode::{Decode, Encode};
use eframe::{self, egui, egui::Context};
use egui_file_dialog::FileDialog;
use file_io::save::{load, StateToSave, DEFAULT_SAVE_FILE};
use gui::navigation::{Page, PageSeq};
use sequence::{seq_from_str, Seq};

use crate::{
    file_io::GenericData,
    gui::{navigation::PageSeqTop, WINDOW_HEIGHT, WINDOW_TITLE, WINDOW_WIDTH},
    pcr::{PcrParams, PolymeraseType},
    primer::TM_TARGET,
    restriction_enzyme::{load_re_library, ReMatch, RestrictionEnzyme},
    sequence::{
        find_orf_matches, seq_to_str, FeatureDirection, FeatureType, ReadingFrame,
        ReadingFrameMatch,
    },
};
use crate::sequence::Feature;

mod features_known;
mod file_io;
mod gui;
mod melting_temp_calcs;
mod pcr;
mod primer;
mod primer_metrics;
mod restriction_enzyme;
mod save_compat;
mod sequence;
mod solution_helper;
mod toxic_proteins;
mod util;

type Color = (u8, u8, u8); // RGB

#[derive(Clone, Copy)]
enum AminoAcid {
    Met,
    Ser,
    Ile,
    Gln,
    His,
    Phe,
    Arg,
    Val,
    Leu,
    Cys,
    Asp,
}

struct PlasmidData {
    // todo : Which seq type? There are many.
    seq_expected: Seq,
    seq_read_assembled: Seq,
    // seq_reads: Vec<SequenceRead>,
}

fn check_seq_integrity(data: &PlasmidData) {}

fn check_toxic_proteins(data: &PlasmidData) {}

fn check_all(data: &PlasmidData) {
    check_seq_integrity(data);
    check_toxic_proteins(data);
}

impl eframe::App for State {
    /// This is the GUI's event loop.
    fn update(&mut self, ctx: &Context, _frame: &mut eframe::Frame) {
        gui::draw(self, ctx);
    }
}

/// Variables for UI fields, for determining PCR parameters.
struct PcrUi {
    pub primer_tm: f32,
    pub product_len: usize,
    pub polymerase_type: PolymeraseType,
    pub num_cycles: u16,
    /// index from primer data. For storing dropdown state.
    pub primer_selected: usize,
}

impl Default for PcrUi {
    fn default() -> Self {
        Self {
            primer_tm: TM_TARGET,
            product_len: 1_000,
            polymerase_type: Default::default(),
            num_cycles: 30,
            primer_selected: 0,
        }
    }
}

#[derive(Clone, Encode, Decode)]
/// Concentrations of common ions in the oglio solution. Affects melting temperature (TM).
/// All values are in milliMolar.
struct IonConcentrations {
    /// Na+ or K+
    pub monovalent: f32,
    /// Mg2+
    pub divalent: f32,
    pub dntp: f32,
    /// Primer concentration, in nM.
    pub primer: f32,
}

impl Default for IonConcentrations {
    fn default() -> Self {
        // todo: Adjust A/R
        Self {
            monovalent: 50.,
            divalent: 1.5,
            dntp: 0.2,
            primer: 25.,
        }
    }
}

#[derive(Default)]
struct StateFeatureAdd {
    // This is in 1-based indexing.
    start_posit: usize,
    end_posit: usize,
    feature_type: FeatureType,
    direction: FeatureDirection,
    label: String,
    color: Option<Color>,
}

/// This Ui struct is used to determine which items on the sequence and map views to show and hide.
struct SeqVisibility {
    /// Show or hide restriction enzymes from the sequence view.
    show_res: bool,
    /// Show and hide primers on
    show_primers: bool,
    /// todo: Show and hide individual features?
    show_features: bool,
    show_reading_frame: bool,
    show_start_stop_codons: bool,
}

impl Default for SeqVisibility {
    fn default() -> Self {
        Self {
            show_res: true,
            show_primers: true,
            show_features: true,
            show_reading_frame: true,
            show_start_stop_codons: false,
        }
    }
}

struct FileDialogs {
    save: FileDialog,
    load: FileDialog,
    import: FileDialog,
    export_fasta: FileDialog,
    export_genbank: FileDialog,
    export_dna: FileDialog,
    selected: Option<PathBuf>,
}

impl Default for FileDialogs {
    fn default() -> Self {
        let save = FileDialog::new()
            // .add_quick_access("Project", |s| {
            //     s.add_path("☆  Examples", "examples");
            // })
            .add_file_filter(
                "PlasCAD files",
                Arc::new(|p| p.extension().unwrap_or_default().to_ascii_lowercase() == "pcad"),
            )
            .default_file_filter("PlasCAD files")
            .id("0");
        // .id("egui_file_dialog");

        let load_ = FileDialog::new()
            .add_file_filter(
                "PlasCAD files",
                Arc::new(|p| p.extension().unwrap_or_default().to_ascii_lowercase() == "pcad"),
            )
            .id("1");

        let import = FileDialog::new()
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
                "FASTA/GB/SnapGene",
                Arc::new(|p| {
                    let ext = p.extension().unwrap_or_default().to_ascii_lowercase();
                    ext == "fasta" || ext == "gb" || ext == "gbk" || ext == "dna"
                }),
            )
            .default_file_filter("FASTA/GB/SnapGene")
            .id("2");

        let export_fasta = FileDialog::new()
            .add_file_filter(
                "FASTA files",
                Arc::new(|p| p.extension().unwrap_or_default().to_ascii_lowercase() == "fasta"),
            )
            .default_file_filter("FASTA files")
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
            .id("4");

        let export_dna = FileDialog::new()
            .add_file_filter(
                "SnapGene DNA files",
                Arc::new(|p| p.extension().unwrap_or_default().to_ascii_lowercase() == "dna"),
            )
            .default_file_filter("SnapGene DNA files")
            .id("5");

        Self {
            save,
            load: load_,
            import,
            export_fasta,
            export_genbank,
            export_dna,
            selected: None,
        }
    }
}

/// Values defined here generally aren't worth saving to file etc.
struct StateUi {
    // todo: Make separate primer cols and primer data; data in state. primer_cols are pre-formatted
    // todo to save computation.
    page: Page,
    page_seq: PageSeq,
    page_seq_top: PageSeqTop,
    seq_insert_input: String,
    seq_vector_input: String,
    seq_input: String,
    pcr: PcrUi,
    feature_add: StateFeatureAdd,
    primer_selected: Option<usize>,
    feature_selected: Option<usize>,
    feature_hover: Option<usize>,
    seq_visibility: SeqVisibility,
    hide_map_feature_editor: bool,
    cursor_pos: Option<(f32, f32)>,
    cursor_seq_i: Option<usize>,
    file_dialogs: FileDialogs,
    /// Show or hide the field to change origin
    show_origin_change: bool,
    new_origin: usize,
}

impl Default for StateUi {
    fn default() -> Self {
        Self {
            page: Default::default(),
            page_seq: Default::default(),
            page_seq_top: Default::default(),
            seq_insert_input: Default::default(),
            seq_vector_input: Default::default(),
            seq_input: Default::default(),
            pcr: Default::default(),
            feature_add: Default::default(),
            primer_selected: None,
            feature_selected: Default::default(),
            feature_hover: Default::default(),
            seq_visibility: Default::default(),
            hide_map_feature_editor: true,
            cursor_pos: None,
            cursor_seq_i: None,
            file_dialogs: Default::default(),
            show_origin_change: false,
            new_origin: 0,
        }
    }
}

pub enum Selection {
    Feature(usize), // index
    Primer(usize),
    None,
}

impl Default for Selection {
    fn default() -> Self {
        Self::None
    }
}

/// Based on GenBank's reference format
#[derive(Default, Clone, Encode, Decode)]
pub struct Reference {
    pub description: String,
    pub authors: Option<String>,
    pub consortium: Option<String>,
    pub title: String,
    pub journal: Option<String>,
    pub pubmed: Option<String>,
    pub remark: Option<String>,
}

/// This struct contains state that does not need to persist between sessesions or saves, but is not
/// a good fit for `StateUi`. This is, generally, calculated data from persistent staet.
#[derive(Default)]
struct StateVolatile {
    restriction_enzyme_sites: Vec<ReMatch>,
    reading_frame_matches: Vec<ReadingFrameMatch>,
}

/// Contains sequence-level metadata.
#[derive(Clone, Default, Encode, Decode)]
pub struct Metadata {
    pub plasmid_name: String,
    pub comments: Vec<String>,
    pub references: Vec<Reference>,
    pub locus: String,
    pub definition: Option<String>,
    pub accession: Option<String>,
    pub version: Option<String>,
    // pub keywords: Vec<String>,
    pub keywords: Option<String>, // todo vec?
    pub source: Option<String>,
    pub organism: Option<String>,
}

/// Note: use of serde traits here and on various sub-structs are for saving and loading.
#[derive(Default)]
struct State {
    ui: StateUi, // Does not need to be saved
    /// Data that is the most fundamental to persistent state, and shared between save formats.
    generic: GenericData,
    // seq: Seq,
    // topology: SeqTopology,
    // features: Vec<Feature>,
    // primers: Vec<Primer>,
    // metadata: Metadata,
    insert_loc: usize,
    // generic: GenericData,// todo: Once your primers are set up.
    ion_concentrations: IonConcentrations,
    pcr: PcrParams,
    restriction_enzyme_lib: Vec<RestrictionEnzyme>, // Does not need to be saved
    selected_item: Selection,
    reading_frame: ReadingFrame,
    volatile: StateVolatile,
}

impl State {
    /// Runs the match serach between primers and sequences. Run this when primers and sequences change.
    pub fn sync_primer_matches(&mut self, primer_i: Option<usize>) {
        let primers = match primer_i {
            Some(i) => &mut self.generic.primers[i..i + 1],
            None => &mut self.generic.primers,
        };

        for primer in primers {
            primer.volatile.matches_seq = primer.match_to_seq(&self.generic.seq);
        }
    }

    pub fn sync_pcr(&mut self) {
        self.pcr = PcrParams::new(&self.ui.pcr);
    }

    /// Identify restriction enzyme sites in the sequence
    pub fn sync_re_sites(&mut self) {
        self.volatile.restriction_enzyme_sites = Vec::new();

        // let seq_comp = seq_complement(seq);

        for (lib_index, re) in self.restriction_enzyme_lib.iter().enumerate() {
            // todo: Use bio lib?
            let seq_len = self.generic.seq.len();
            for i in 0..seq_len {
                if i + re.seq.len() + 1 >= seq_len {
                    continue;
                }

                if re.seq == self.generic.seq[i..i + re.seq.len()] {
                    self.volatile.restriction_enzyme_sites.push(ReMatch {
                        lib_index,
                        seq_index: i,
                        // direction: PrimerDirection::Forward,
                    });
                }
                // todo: Simpler way?
                // todo: Evaluate if this is what you want.
                // if re.seq == seq_comp[i..i + re.seq.len()] {
                //     self.restriction_enzyme_sites.push(ReMatch {
                //         lib_index: i_re,
                //         seq_index: i,
                //         direction: PrimerDirection::Reverse,
                //     });
                // }
            }
        }

        // This sorting aids in our up/down label alternation in the display.
        self.volatile
            .restriction_enzyme_sites
            .sort_by(|a, b| a.seq_index.cmp(&b.seq_index));
    }

    pub fn sync_reading_frame(&mut self) {
        self.volatile.reading_frame_matches =
            find_orf_matches(&self.generic.seq, self.reading_frame);
    }

    pub fn sync_primer_metrics(&mut self) {
        for primer in &mut self.generic.primers {
            primer.volatile.sequence_input = seq_to_str(&primer.sequence);

            if let Some(metrics) = &mut primer.volatile.metrics {
                metrics.update_scores();
            }
            if primer.volatile.metrics.is_none() {
                primer.run_calcs(&self.ion_concentrations);
            }
        }
    }

    /// Update the combined SLIC vector + insert sequence.
    pub fn sync_cloning_product(&mut self) {
        let seq_vector = seq_from_str(&self.ui.seq_vector_input);
        let seq_insert = seq_from_str(&self.ui.seq_insert_input);

        if self.insert_loc + 1 > seq_vector.len() {
            eprintln!(
                "Error creating cloning insert: insert loc {} is greater than vector len {}",
                self.insert_loc,
                seq_vector.len()
            );
            return;
        }

        self.generic.seq.clone_from(&seq_vector);
        self.generic
            .seq
            .splice(self.insert_loc..self.insert_loc, seq_insert.iter().cloned());

        self.sync_primer_matches(None);
    }

    /// Run this when the sequence changes.
    pub fn sync_seq_related(&mut self, primer_i: Option<usize>) {
        self.sync_primer_matches(primer_i);
        self.sync_re_sites();
        self.sync_reading_frame();

        self.ui.seq_input = seq_to_str(&self.generic.seq);
    }

    /// Load state from a (our format) file.
    pub fn load(path: &str) -> Self {
        let state_loaded: io::Result<StateToSave> = load(path);
        let mut result = match state_loaded {
            Ok(s) => s.to_state(),
            Err(_) => Default::default(),
        };

        result.restriction_enzyme_lib = load_re_library();

        result.sync_pcr();
        result.sync_primer_metrics();
        result.sync_seq_related(None);
        result.ui.seq_input = seq_to_str(&result.generic.seq);

        result
    }
}

fn main() {
    let state = State::load(DEFAULT_SAVE_FILE);

    let icon_bytes: &[u8] = include_bytes!("resources/icon.png");
    let icon_data = eframe::icon_data::from_png_bytes(icon_bytes);

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([WINDOW_WIDTH, WINDOW_HEIGHT])
            .with_icon(icon_data.unwrap()),
        follow_system_theme: false,
        ..Default::default()
    };

    eframe::run_native(WINDOW_TITLE, options, Box::new(|_cc| Ok(Box::new(state)))).unwrap();
}
