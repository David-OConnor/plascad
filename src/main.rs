// Disables the terminal window. Use this for releases, but disable when debugging.
// #![windows_subsystem = "windows"]

use std::{
    io,
    path::{Path, PathBuf},
    sync::Arc,
};

use bincode::{config, Decode, Encode};
// use bio::{
//     bio_types::sequence::{Sequence, SequenceRead},
//     data_structures::fmindex::FMIndexable,
//     io::fastq::FastqRead,
// };
use eframe::{self, egui, egui::Context};
// use egui_file::FileDialog;
use egui_file_dialog::FileDialog;
use gui::navigation::{Page, PageSeq};
use primer::PrimerData;
use save::load;
use sequence::{seq_from_str, Seq};

use crate::{
    gui::{navigation::PageSeqTop, WINDOW_HEIGHT, WINDOW_TITLE, WINDOW_WIDTH},
    pcr::{PcrParams, PolymeraseType},
    primer::{PrimerDirection, TM_TARGET},
    restriction_enzyme::{load_re_library, NucleotideGeneral::N, ReMatch, RestrictionEnzyme},
    save::{StateToSave, DEFAULT_SAVE_FILE},
    sequence::{
        seq_to_str, Feature, FeatureDirection, FeatureType, Nucleotide,
        Nucleotide::{A, G, T},
        ReadingFrame, ReadingFrameMatch, SeqTopology,
    },
};

mod features_known;
mod genbank_parse;
mod gui;
mod melting_temp_calcs;
mod pcr;
mod primer;
mod primer_metrics;
mod restriction_enzyme;
mod save;
mod save_compat;
mod sequence;
mod snapgene_parse;
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
    // page_primer: PagePrimer,
    page_seq: PageSeq,
    page_seq_top: PageSeqTop,
    seq_insert_input: String,
    seq_vector_input: String,
    seq_input: String,
    pcr: PcrUi,
    feature_add: StateFeatureAdd,
    // pcr_primer: Option<usize>, // primer index, if primer count > 0.
    // pcr_primer: usize, // primer index
    primer_selected: Option<usize>,
    // todo: suitstruct for show/hide A/R
    // /// Hide the primer addition/creation/QC panel.
    // hide_primer_table: bool,
    seq_visibility: SeqVisibility,
    hide_map_feature_editor: bool,
    cursor_pos: Option<(f32, f32)>,
    cursor_seq_i: Option<usize>,
    file_dialogs: FileDialogs,
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
            // pcr_primer: 0,
            primer_selected: None,
            // hide_primer_table: false,
            seq_visibility: Default::default(),
            hide_map_feature_editor: true,
            cursor_pos: None,
            cursor_seq_i: None,
            file_dialogs: Default::default(),
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

/// This struct contains state that does not need to persist between sessesions or saves, but is not
/// a good fit for `StateUi`. This is, generally, calculated data from persistent staet.
#[derive(Default)]
struct StateVolatile {
    restriction_enzyme_sites: Vec<ReMatch>,
    reading_frame_matches: Vec<ReadingFrameMatch>,
}

/// Note: use of serde traits here and on various sub-structs are for saving and loading.
#[derive(Default)]
struct State {
    ui: StateUi, // Does not need to be saved
    primer_data: Vec<PrimerData>,
    // /// Insert and vector are for SLIC and FC.
    // seq_insert: Seq,
    // seq_vector: Seq,
    seq: Seq,
    // /// Amplicon is for basic PCR.
    // seq_amplicon: Seq,
    insert_loc: usize,
    // /// These limits for choosing the insert location may be defined by the vector's promoter, RBS etc.
    // insert_location_5p_limit: usize,
    // insert_location_3p_limit: usize,
    ion_concentrations: IonConcentrations,
    pcr: PcrParams,
    restriction_enzyme_lib: Vec<RestrictionEnzyme>, // Does not need to be saved
    features: Vec<Feature>,
    plasmid_name: String,
    topology: SeqTopology,
    selected_item: Selection,
    reading_frame: ReadingFrame,
    volatile: StateVolatile,
}

impl State {
    /// Runs the match serach between primers and sequences. Run this when primers and sequences change.
    pub fn sync_primer_matches(&mut self, primer_i: Option<usize>) {
        let p_list = match primer_i {
            Some(i) => &mut self.primer_data[i..i + 1],
            None => &mut self.primer_data,
        };

        for p_data in p_list {
            p_data.matches_seq = p_data.primer.match_to_seq(&self.seq);
            // p_data.matches_insert = p_data.primer.match_to_seq(&self.seq_insert);
            // p_data.matches_vector = p_data.primer.match_to_seq(&self.seq_vector);
            // p_data.matches_vector_with_insert =
            //     p_data.primer.match_to_seq(&self.seq);
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
            for i in 0..self.seq.len() {
                if i + re.seq.len() + 1 >= self.seq.len() {
                    continue;
                }

                if re.seq == self.seq[i..i + re.seq.len()] {
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
        const START_CODON: [Nucleotide; 3] = [A, T, G];

        self.volatile.reading_frame_matches = Vec::new();

        let seq = match self.reading_frame {
            ReadingFrame::Fwd0 | ReadingFrame::Fwd1 | ReadingFrame::Fwd2 => (),
            _ => (),
        };
    }

    pub fn sync_primer_metrics(&mut self) {
        for primer in &mut self.primer_data {
            if let Some(metrics) = &mut primer.metrics {
                metrics.update_scores();
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

        self.seq = seq_vector.clone(); // Clone the original vector
        self.seq
            .splice(self.insert_loc..self.insert_loc, seq_insert.iter().cloned());

        self.sync_primer_matches(None);
    }

    /// Run this when the sequence changes.
    pub fn sync_seq_related(&mut self, primer_i: Option<usize>) {
        self.sync_primer_matches(primer_i);
        self.sync_re_sites();
        self.sync_reading_frame();
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
        result.ui.seq_input = seq_to_str(&result.seq);

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
