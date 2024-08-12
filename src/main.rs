// Disables the terminal window on Windows, in release mode.
#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

// todo: Build a database of feature sequences. You can find GenBank etc files online (addGene, and other sources)
// todo: and parse common features.

// todo: Break out Generic into its own mod?

use std::{
    io,
    path::{Path, PathBuf},
    str::FromStr,
    sync::Arc,
};

use bincode::{Decode, Encode};
use eframe::{self, egui, egui::Context};
use egui_file_dialog::{FileDialog, FileDialogConfig};
use file_io::save::{load, StateToSave, DEFAULT_SAVE_FILE};
use gui::navigation::{Page, PageSeq};
use primer::IonConcentrations;
use sequence::{seq_from_str, Seq};

use crate::{
    file_io::{
        save::{
            StateUiToSave, DEFAULT_DNA_FILE, DEFAULT_FASTA_FILE, DEFAULT_GENBANK_FILE,
            DEFAULT_PREFS_FILE,
        },
        GenericData,
    },
    gui::{navigation::PageSeqTop, WINDOW_HEIGHT, WINDOW_TITLE, WINDOW_WIDTH},
    pcr::{PcrParams, PolymeraseType},
    primer::TM_TARGET,
    restriction_enzyme::{load_re_library, ReMatch, RestrictionEnzyme},
    sequence::{
        find_orf_matches, seq_to_str, Feature, FeatureDirection, FeatureType, ReadingFrame,
        ReadingFrameMatch,
    },
};

mod feature_db_load;
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
    fn update(&mut self, ctx: &Context, frame: &mut eframe::Frame) {
        gui::draw(self, ctx);
    }
}

#[derive(Clone, Encode, Decode)]
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

#[derive(Default, Encode, Decode)]
struct StateFeatureAdd {
    // This is in 1-based indexing.
    start_posit: usize,
    end_posit: usize,
    feature_type: FeatureType,
    direction: FeatureDirection,
    label: String,
    color: Option<Color>,
}

#[derive(Clone, Encode, Decode)]
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
    // load: FileDialog,
    load: FileDialog,
    export_fasta: FileDialog,
    export_genbank: FileDialog,
    export_dna: FileDialog,
    cloning_load: FileDialog,
    // todo: What do we use this for?
    // selected: Option<PathBuf>,
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
            "PCAD/FASTA/GB/SG",
            Arc::new(|p| {
                let ext = p.extension().unwrap_or_default().to_ascii_lowercase();
                ext == "pcad" || ext == "fasta" || ext == "gb" || ext == "gbk" || ext == "dna"
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
            .default_file_name(DEFAULT_SAVE_FILE)
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

#[derive(Default)]
struct CloningInsertData {
    /// We use this list to store a list of features to clone an insert from, loaded from a file.
    pub features_loaded: Vec<Feature>,
    /// `cloning_ins_features_loaded` indexes reference this sequence.
    pub seq_loaded: Seq,
    pub feature_selected: Option<usize>,
    pub seq_insert: Seq,
    pub seq_input: String,
}

/// Values defined here generally aren't worth saving to file etc.
struct StateUi {
    // todo: Make separate primer cols and primer data; data in state. primer_cols are pre-formatted
    // todo to save computation.
    page: Page,
    page_seq: PageSeq,
    page_seq_top: PageSeqTop,
    seq_input: String,
    pcr: PcrUi,
    feature_add: StateFeatureAdd,
    primer_selected: Option<usize>, // primer page only.
    feature_hover: Option<usize>,
    selected_item: Selection,
    seq_visibility: SeqVisibility,
    hide_map_feature_editor: bool,
    /// Mouse cursor
    cursor_pos: Option<(f32, f32)>,
    /// Mouse cursor
    cursor_seq_i: Option<usize>,
    file_dialogs: FileDialogs,
    /// Show or hide the field to change origin
    show_origin_change: bool,
    new_origin: usize,
    /// Text-editing cursor. Used for editing on the sequence view. Chars typed
    /// will be inserted after this index. This index is 0-based.
    text_cursor_i: Option<usize>,
    /// We store if we've clicked somewhere separately from the action, as getting a sequence index
    /// from cursor positions may be decoupled, and depends on the view.
    click_pending_handle: bool,
    cloning_insert: CloningInsertData,
}

impl Default for StateUi {
    fn default() -> Self {
        Self {
            page: Default::default(),
            page_seq: Default::default(),
            page_seq_top: Default::default(),
            seq_input: Default::default(),
            pcr: Default::default(),
            feature_add: Default::default(),
            primer_selected: None,
            feature_hover: Default::default(),
            selected_item: Selection::None,
            seq_visibility: Default::default(),
            hide_map_feature_editor: true,
            cursor_pos: None,
            cursor_seq_i: None,
            file_dialogs: Default::default(),
            show_origin_change: false,
            new_origin: 0,
            text_cursor_i: None,
            click_pending_handle: false,
            cloning_insert: Default::default(),
        }
    }
}

#[derive(Clone, Copy, Encode, Decode)]
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
    reading_frame: ReadingFrame,
    volatile: StateVolatile,
    /// Used to determine how the save function works, among other things.
    path_loaded: Option<PathBuf>,
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

    /// Upddate this sequence by inserting a sequence of interest; the new sequence is the cloning product.
    pub fn sync_cloning_product(&mut self) {
        let seq_vector = &mut self.generic.seq;
        let seq_insert = &self.ui.cloning_insert.seq_insert;

        if self.insert_loc + 1 > seq_vector.len() {
            eprintln!(
                "Error creating cloning insert: insert loc {} is greater than vector len {}",
                self.insert_loc,
                seq_vector.len()
            );
            return;
        }

        // self.generic.seq.clone_from(&seq_vector);
        seq_vector.splice(self.insert_loc..self.insert_loc, seq_insert.iter().cloned());

        // todo: YOu may run into off-by-one issues on insert loc here; adjust A/R.
        // Now, you have to update features affected by this insertion, shifting them right A/R.
        for feature in &mut self.generic.features {
            // Handle the case where the insert occurs over a feature. Do this before shifting features.
            if feature.index_range.0 < self.insert_loc && feature.index_range.1 > self.insert_loc {
                // todo: Divide into two features? For now, we are just trimming.
                feature.index_range.1 = self.insert_loc - 1;
            }

            if feature.index_range.0 > self.insert_loc {
                feature.index_range.0 += seq_insert.len();
            }

            if feature.index_range.1 > self.insert_loc {
                feature.index_range.1 += seq_insert.len();
            }
        }

        self.sync_seq_related(None);
    }

    /// Run this when the sequence changes.
    pub fn sync_seq_related(&mut self, primer_i: Option<usize>) {
        self.sync_primer_matches(primer_i);
        self.sync_re_sites();
        self.sync_reading_frame();

        self.ui.seq_input = seq_to_str(&self.generic.seq);
    }

    /// Load state from a (our format) file.
    pub fn load(path: &Path, prefs_path: &Path) -> Self {
        let state_loaded: io::Result<StateToSave> = load(path);
        let mut result = match state_loaded {
            Ok(s) => s.to_state(),
            Err(_) => Default::default(),
        };

        let ui_loaded: io::Result<StateUiToSave> = load(prefs_path);

        result.ui = match ui_loaded {
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
    let state = State::load(
        &PathBuf::from_str(DEFAULT_SAVE_FILE).unwrap(),
        &PathBuf::from_str(DEFAULT_PREFS_FILE).unwrap(),
    );

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
