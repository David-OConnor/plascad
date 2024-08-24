// Disables the terminal window on Windows, in release mode.
#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

// todo: Build a database of feature sequences. You can find GenBank etc files online (addGene, and other sources)
// todo: and parse common features.

// todo: Break out Generic into its own mod?

// Reading frame: Guess the frame, and truncate the start based on CodingRegion and Gene feature types?
use std::{
    env, io,
    path::{Path, PathBuf},
    str::FromStr,
    sync::Arc,
};

use bincode::{Decode, Encode};
use cloning::CloningInsertData;
use copypasta::{ClipboardContext, ClipboardProvider};
use eframe::{self, egui, egui::Context};
use egui_file_dialog::{FileDialog, FileDialogConfig};
use file_io::save::{load, StateToSave, DEFAULT_SAVE_FILE};
use gui::navigation::{Page, PageSeq};
use primer::IonConcentrations;
use protein::Protein;
use sequence::Seq;

use crate::{
    amino_acids::AaIdent,
    file_io::{
        save::{
            save, StateUiToSave, DEFAULT_DNA_FILE, DEFAULT_FASTA_FILE, DEFAULT_GENBANK_FILE,
            DEFAULT_PREFS_FILE,
        },
        GenericData,
    },
    gui::{navigation::PageSeqTop, save::load_import, WINDOW_HEIGHT, WINDOW_TITLE, WINDOW_WIDTH},
    pcr::{PcrParams, PolymeraseType},
    portions::PortionsState,
    primer::TM_TARGET,
    protein::sync_cr_orf_matches,
    restriction_enzyme::{find_re_matches, load_re_library, ReMatch, RestrictionEnzyme},
    sequence::{
        find_orf_matches, find_search_matches, seq_to_str, FeatureDirection, FeatureType,
        Nucleotide, ReadingFrame, ReadingFrameMatch, SearchMatch, MIN_SEARCH_LEN,
    },
    tags::TagMatch,
    util::RangeIncl,
};

mod amino_acids;
mod cloning;
mod feature_db_load;
mod file_io;
mod gui;
mod melting_temp_calcs;
mod pcr;
mod portions;
mod primer;
mod primer_metrics;
mod protein;
mod restriction_enzyme;
mod save_compat;
mod sequence;
mod solution_helper;
mod tags;
mod toxic_proteins;
mod util;

type Color = (u8, u8, u8); // RGB

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
    /// Index from primer data for the load-from-primer system. For storing dropdown state.
    pub primer_selected: usize,
    /// These are for the PCR product generation
    pub primer_fwd: usize,
    pub primer_rev: usize,
}

impl Default for PcrUi {
    fn default() -> Self {
        Self {
            primer_tm: TM_TARGET,
            product_len: 1_000,
            polymerase_type: Default::default(),
            num_cycles: 30,
            primer_selected: 0,
            primer_fwd: 0,
            primer_rev: 0,
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
}

impl Default for SeqVisibility {
    fn default() -> Self {
        Self {
            show_res: true,
            show_primers: true,
            show_features: true,
            show_reading_frame: false,
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
    // primer_selected: Option<usize>, // primer page only.
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
    /// We use this for selecting features from the seq view
    dblclick_pending_handle: bool,
    cloning_insert: CloningInsertData,
    /// Volatile; computed dynamically based on window size.
    nt_chars_per_row: usize,
    search_input: String,
    /// Used to trigger a search focus on hitting ctrl+f
    highlight_search_input: bool,
    /// Activated when the user selects the search box; disables character insertion.
    text_edit_active: bool,
    /// This is used for selecting nucleotides on the sequence viewer.
    dragging: bool,
    /// 1-based indexing.
    text_selection: Option<RangeIncl>,
    quick_feature_add_name: String,
    quick_feature_add_dir: FeatureDirection,
    aa_ident_disp: AaIdent,
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
            feature_hover: Default::default(),
            selected_item: Default::default(),
            seq_visibility: Default::default(),
            hide_map_feature_editor: true,
            cursor_pos: Default::default(),
            cursor_seq_i: Default::default(),
            file_dialogs: Default::default(),
            show_origin_change: Default::default(),
            new_origin: Default::default(),
            text_cursor_i: Some(0),
            click_pending_handle: Default::default(),
            dblclick_pending_handle: Default::default(),
            cloning_insert: Default::default(),
            nt_chars_per_row: Default::default(),
            search_input: Default::default(),
            highlight_search_input: Default::default(),
            text_edit_active: Default::default(),
            dragging: Default::default(),
            text_selection: Default::default(),
            quick_feature_add_name: Default::default(),
            quick_feature_add_dir: Default::default(),
            aa_ident_disp: AaIdent::ThreeLetters,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Encode, Decode)]
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
    tag_matches: Vec<TagMatch>,
    search_matches: Vec<SearchMatch>,
    /// Used for automatically determining which reading frame to use, and the full frame,
    /// for a given coding-region feature.
    cr_orf_matches: Vec<(usize, ReadingFrameMatch)>,
}

/// Note: use of serde traits here and on various sub-structs are for saving and loading.
#[derive(Default)]
struct State {
    ui: StateUi,
    /// Data that is the most fundamental to persistent state, and shared between save formats.
    generic: GenericData,
    cloning_insert_loc: usize,
    ion_concentrations: IonConcentrations,
    pcr: PcrParams,
    restriction_enzyme_lib: Vec<RestrictionEnzyme>, // Does not need to be saved
    reading_frame: ReadingFrame,
    volatile: StateVolatile,
    /// Used to determine how the save function works, among other things.
    path_loaded: Option<PathBuf>,
    search_seq: Seq,
    portions: PortionsState,
    proteins: Vec<Protein>,
}

impl State {
    /// Reset data; we currently use this for making "new" data.
    pub fn reset(&mut self) {
        self.generic = Default::default();
        self.volatile = Default::default();

        self.ui.cursor_pos = None;
        self.ui.cursor_seq_i = None;
        self.ui.text_cursor_i = Some(0); // todo: For now; having trouble with cursor on empty seq
        self.ui.seq_input = String::new();
    }

    /// Load UI and related data not related to a specific sequence.
    pub fn load_prefs(&mut self, path: &Path) {
        let ui_loaded: io::Result<StateUiToSave> = load(path);

        if let Ok(ui) = ui_loaded {
            self.ui = ui.to_state();
        }
    }

    pub fn save_prefs(&self) {
        if let Err(e) = save(
            &PathBuf::from(DEFAULT_PREFS_FILE),
            &StateUiToSave::from_state(&self.ui),
        ) {
            eprintln!("Error saving prefs: {e}");
        }
    }

    /// Runs the match serach between primers and sequences. Run this when primers and sequences change.
    pub fn sync_primer_matches(&mut self, primer_i: Option<usize>) {
        let primers = match primer_i {
            Some(i) => &mut self.generic.primers[i..=i],
            None => &mut self.generic.primers,
        };

        for primer in primers {
            primer.volatile.matches = primer.match_to_seq(&self.generic.seq);
        }
    }

    pub fn sync_pcr(&mut self) {
        self.pcr = PcrParams::new(&self.ui.pcr);
    }

    /// Identify restriction enzyme sites in the sequence
    pub fn sync_re_sites(&mut self) {
        self.volatile.restriction_enzyme_sites = Vec::new();

        self.volatile
            .restriction_enzyme_sites
            .append(&mut find_re_matches(
                &self.generic.seq,
                &self.restriction_enzyme_lib,
            ));

        // This sorting aids in our up/down label alternation in the display.
        self.volatile
            .restriction_enzyme_sites
            .sort_by(|a, b| a.seq_index.cmp(&b.seq_index));
    }

    pub fn sync_reading_frame(&mut self) {
        self.volatile.reading_frame_matches =
            find_orf_matches(&self.generic.seq, self.reading_frame);
    }

    pub fn sync_search(&mut self) {
        if self.search_seq.len() >= MIN_SEARCH_LEN {
            self.volatile.search_matches = find_search_matches(&self.generic.seq, &self.search_seq);
        } else {
            self.volatile.search_matches = Vec::new();
        }
    }

    pub fn sync_portions(&mut self) {
        for sol in &mut self.portions.solutions {
            sol.calc_amounts();
        }
    }

    pub fn sync_primer_metrics(&mut self) {
        for primer in &mut self.generic.primers {
            primer.run_calcs(&self.ion_concentrations);
            //
            // primer.volatile.sequence_input = seq_to_str(&primer.sequence);
            //
            //
            // if let Some(metrics) = &mut primer.volatile.metrics {
            //     let dual_ended = match primer.volatile.tune_setting {
            //         TuneSetting::Both(_) => true,
            //         _ => false,
            //     };
            //
            //     metrics.update_scores(dual_ended);
            // }
            // if primer.volatile.metrics.is_none() {
            //     primer.run_calcs(&self.ion_concentrations);
            // }
        }
    }

    /// Upddate this sequence by inserting a sequence of interest. `insert_loc` uses 0-based indexing.
    pub fn insert_nucleotides(&mut self, insert: &[Nucleotide], insert_loc: usize) {
        let seq_vector = &mut self.generic.seq;

        let insert_i = insert_loc - 1; // 1-based indexing.

        if insert_i > seq_vector.len() {
            eprintln!(
                "Error inserting nucleotides: insert loc {} is greater than vector len {}",
                insert_loc,
                seq_vector.len()
            );
            return;
        }

        seq_vector.splice(insert_i..insert_i, insert.iter().cloned());

        // Now, you have to update features affected by this insertion, shifting them right A/R.
        for feature in &mut self.generic.features {
            if feature.range.start > insert_i {
                feature.range.start = feature.range.start + insert.len();
            }

            if feature.range.end > insert_i {
                feature.range.end = feature.range.end + insert.len();
            }
        }

        self.sync_seq_related(None);
    }

    /// One-based indexing. Similar to `insert_nucleotides`.
    pub fn remove_nucleotides(&mut self, range: RangeIncl) {
        let seq = &mut self.generic.seq;
        if range.end + 1 > seq.len() {
            return;
        }

        let count = range.len();

        seq.drain(range.start..=range.end);

        // Now, you have to update features affected by this insertion, shifting them left A/R.
        for feature in &mut self.generic.features {
            if feature.range.start > range.end {
                feature.range.start = feature.range.start - count;
            }
            if feature.range.end > range.end {
                feature.range.end = feature.range.end - count;
            }
        }

        self.sync_seq_related(None);
    }

    /// Run this when the sequence changes.
    pub fn sync_seq_related(&mut self, primer_i: Option<usize>) {
        self.sync_primer_matches(primer_i);
        self.sync_re_sites();
        self.sync_reading_frame();
        self.sync_search();

        sync_cr_orf_matches(self);

        self.ui.seq_input = seq_to_str(&self.generic.seq);
    }

    /// Load state from a (our format) file.
    pub fn load(path: &Path, prefs_path: &Path) -> Self {
        let state_loaded: io::Result<StateToSave> = load(path);
        let mut result = match state_loaded {
            Ok(s) => s.to_state(),
            Err(_) => Default::default(),
        };

        result.load_prefs(prefs_path);

        result.restriction_enzyme_lib = load_re_library();

        result.sync_pcr();
        result.sync_primer_metrics();
        result.sync_seq_related(None);
        result.ui.seq_input = seq_to_str(&result.generic.seq);
        result.sync_portions();

        result
    }

    /// Copy the sequence of the selected text selection, feature or primer to the clipboard, if applicable.
    pub fn copy_seq(&self) {
        // Text selection takes priority.
        if let Some(selection) = &self.ui.text_selection {
            if let Some(seq) = selection.index_seq(&self.generic.seq) {
                let mut ctx = ClipboardContext::new().unwrap();
                ctx.set_contents(seq_to_str(seq)).unwrap();
            }
            return;
        }

        match self.ui.selected_item {
            Selection::Feature(i) => {
                let feature = &self.generic.features[i];
                // todo: Why -1 not working?
                if let Some(seq) = feature.range.index_seq(&self.generic.seq) {
                    let mut ctx = ClipboardContext::new().unwrap();
                    ctx.set_contents(seq_to_str(seq)).unwrap();
                }
            }
            Selection::Primer(i) => {
                let primer = &self.generic.primers[i];

                let mut ctx = ClipboardContext::new().unwrap();
                ctx.set_contents(seq_to_str(&primer.sequence)).unwrap();
            }
            _ => (),
        }
    }
}

fn main() {
    let mut window_title_initial = WINDOW_TITLE.to_owned();
    let path = {
        let mut r = PathBuf::from_str(DEFAULT_SAVE_FILE).unwrap();

        // Windows and possibly other operating systems, if attempting to use your program to natively
        // open a file type, will use command line argumentse to indicate this. Determine if the program
        // is being launched this way, and if so, open the file.
        let args: Vec<String> = env::args().collect();
        if args.len() > 1 {
            let temp = &args[1];
            r = PathBuf::from_str(&temp).unwrap();

            // Just the filename and extension.
            window_title_initial = r
                .file_name()
                .and_then(|name| name.to_str())
                .map(|name_str| name_str.to_string())
                .unwrap();
        }

        r
    };

    let mut state = State::default();
    load_import(&mut state, &path);

    state.load_prefs(&PathBuf::from_str(DEFAULT_PREFS_FILE).unwrap());

    if window_title_initial != WINDOW_TITLE {
        state.path_loaded = Some(path);
    }

    let icon_bytes: &[u8] = include_bytes!("resources/icon.png");
    let icon_data = eframe::icon_data::from_png_bytes(icon_bytes);

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([WINDOW_WIDTH, WINDOW_HEIGHT])
            .with_icon(icon_data.unwrap()),
        follow_system_theme: false,
        ..Default::default()
    };

    eframe::run_native(
        &window_title_initial,
        options,
        Box::new(|_cc| Ok(Box::new(state))),
    )
    .unwrap();
}
