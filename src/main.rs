// Disables the terminal window on Windows, in release mode.
#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

// todo: Build a database of feature sequences. You can find GenBank etc files online (addGene, and other sources)
// todo: and parse common features.

// todo: Break out Generic into its own mod?

// todo:
// Focus on tools that allow you to conveniently design plasmid seqs based on source vectors, REs etc. It will make primers,
// choose how to combine sequences etc. so, more towards realistic, product-first workflows.
//
// This will make more sense when you are designing new products.
// For example: consider a lib of Addgene's generic vectors for E. Coli.
//
// The input: your target product: Output: as much we can automate as possible.

// Reading frame: Guess the frame, and truncate the start based on CodingRegion and Gene feature types?
use std::{
    env, io,
    path::{Path, PathBuf},
    str::FromStr,
    sync::Arc,
    time::Instant,
};

use bincode::{Decode, Encode};
use cloning::{AutocloneStatus, CloningInsertData};
use copypasta::{ClipboardContext, ClipboardProvider};
use eframe::{self, egui, egui::Context};
use egui_file_dialog::{FileDialog, FileDialogConfig};
use file_io::save::{load, load_import, StateToSave, QUICKSAVE_FILE};
use gui::navigation::{Page, PageSeq};
use primer::IonConcentrations;
use protein::Protein;
use reading_frame::{find_orf_matches, ReadingFrame, ReadingFrameMatch};
use sequence::Seq;

use crate::{
    amino_acids::AaIdent,
    backbones::{load_backbone_library, Backbone, BackboneFilters},
    cloning::BackboneSelected,
    file_io::{
        save::{
            save, StateUiToSave, DEFAULT_DNA_FILE, DEFAULT_FASTA_FILE, DEFAULT_GENBANK_FILE,
            DEFAULT_PREFS_FILE,
        },
        GenericData,
    },
    gui::{navigation::PageSeqTop, WINDOW_HEIGHT, WINDOW_WIDTH},
    ligation::LigationFragment,
    pcr::{PcrParams, PolymeraseType},
    portions::PortionsState,
    primer::TM_TARGET,
    protein::{proteins_from_seq, sync_cr_orf_matches},
    restriction_enzyme::{find_re_matches, load_re_library, ReMatch, RestrictionEnzyme},
    sequence::{
        find_search_matches, seq_to_str, FeatureDirection, FeatureType, Nucleotide, SearchMatch,
        MIN_SEARCH_LEN,
    },
    tags::TagMatch,
    util::{get_window_title, RangeIncl},
};

mod amino_acids;
mod backbones;
mod cloning;
mod external_websites;
mod feature_db_load;
mod file_io;
mod gui;
mod ligation;
mod melting_temp_calcs;
mod pcr;
mod portions;
mod primer;
mod primer_metrics;
mod protein;
mod reading_frame;
mod restriction_enzyme;
mod save_compat;
mod sequence;
mod solution_helper;
mod tags;
mod toxic_proteins;
mod util;

type Color = (u8, u8, u8); // RGB

// todo: Eventually, implement a system that automatically checks for changes, and don't
// todo save to disk if there are no changes.
const PREFS_SAVE_INTERVAL: u64 = 60; // Save user preferences this often, in seconds.

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
    /// This is the GUI's event loop. This also handles periodically saving preferences to disk.
    /// Note that preferences are only saved if the window is active, ie mouse movement in it or similar.
    fn update(&mut self, ctx: &Context, _frame: &mut eframe::Frame) {
        // Note that if the window is
        static mut LAST_PREF_SAVE: Option<Instant> = None;

        let now = Instant::now();

        unsafe {
            if let Some(last_save) = LAST_PREF_SAVE {
                if (now - last_save).as_secs() > PREFS_SAVE_INTERVAL {
                    LAST_PREF_SAVE = Some(now);
                    self.save_prefs()
                }
            } else {
                // Initialize LAST_PREF_SAVE the first time it's accessed
                LAST_PREF_SAVE = Some(now);
            }
        }

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
            .default_file_name(QUICKSAVE_FILE)
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

/// UI state for restriction enzymes.
struct ReUi {
    /// Inner: RE name
    /// todo: This is a trap for multiple tabs.
    // res_selected: Vec<String>,
    res_selected: Vec<RestrictionEnzyme>,
    /// Which tabs' sequences to digest. Note that the index of this vec doesn't matter; the values, which
    /// point to indices elsewhere, does.
    tabs_selected: Vec<usize>,
    unique_cutters_only: bool,
    /// No blunt ends; must produce overhangs.
    sticky_ends_only: bool,
    /// Only show REs that are present in at least two sequences.
    multiple_seqs: bool,
}

impl Default for ReUi {
    fn default() -> Self {
        Self {
            res_selected: Default::default(),
            tabs_selected: Default::default(),
            unique_cutters_only: true,
            sticky_ends_only: false,
            multiple_seqs: true,
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
    seq_input: String, // todo: Consider moving this to volatile.
    pcr: PcrUi,
    feature_add: StateFeatureAdd,
    // primer_selected: Option<usize>, // primer page only.
    feature_hover: Option<usize>, // todo: Apply similar enum logic to selection: Allow Primer::, Feature::, or None::
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
    // todo: Protein ui A/R
    aa_ident_disp: AaIdent,
    pdb_error_received: bool,
    re: ReUi,
    backbone_filters: BackboneFilters,
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
            pdb_error_received: false,
            re: Default::default(),
            backbone_filters: Default::default(),
        }
    }
}

#[derive(Clone, Copy, PartialEq, Debug, Encode, Decode)]
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
    restriction_enzyme_matches: Vec<ReMatch>,
    re_digestion_products: Vec<LigationFragment>,
    reading_frame_matches: Vec<ReadingFrameMatch>,
    tag_matches: Vec<TagMatch>,
    search_matches: Vec<SearchMatch>,
    /// Used for automatically determining which reading frame to use, and the full frame,
    /// for a given coding-region feature.
    cr_orf_matches: Vec<(usize, ReadingFrameMatch)>,
    proteins: Vec<Protein>,
}

struct CloningState {
    backbone_selected: BackboneSelected,
    /// Note: This is only used currently if using the opened file as the BB; otherwise
    /// we use a ref to the library.
    backbone: Option<Backbone>, // todo: Other options: Store a ref; use in bb_selected instead of an index.
    // res_matched:Vec<usize>, // todo: A/R
    res_common: Vec<RestrictionEnzyme>,
    re_matches_vec_common: Vec<ReMatch>,
    re_matches_insert_common: Vec<ReMatch>,
    status: AutocloneStatus, // todo: Should this be an option?
    insert_loc: usize,
    // /// We use this, for example, for displaying a linear map based on a library backbone.
    // backbone_data: Option<GenericData>,
}

impl Default for CloningState {
    fn default() -> Self {
        Self {
            insert_loc: 1,
            backbone_selected: Default::default(),
            backbone: Default::default(),
            res_common: Default::default(),
            re_matches_vec_common: Default::default(),
            re_matches_insert_common: Default::default(),
            status: Default::default(),
            // backbone_data: Default::default(),
        }
    }
}

/// Note: use of serde traits here and on various sub-structs are for saving and loading.
// #[derive(Default)]
struct State {
    ui: StateUi,
    /// Ie tab, file etc.
    active: usize,
    // todo: Consider grouping generic, path_loaded, portions, and similar in a single vec.
    // todo: Do that after your initial tab approach works.
    /// Data that is the most fundamental to persistent state, and shared between save formats.
    /// The index corresponds to `active`.`
    generic: Vec<GenericData>,
    /// Used to determine how the save function works, among other things.
    /// Index corresponds to `active`.
    path_loaded: Vec<Option<PathBuf>>,
    /// Index corresponds to `active`.
    portions: Vec<PortionsState>,
    /// Index corresponds to `active`.
    volatile: Vec<StateVolatile>,
    /// Used for PCR. Index corresponds to `active`.
    ion_concentrations: Vec<IonConcentrations>,
    pcr: PcrParams,
    restriction_enzyme_lib: Vec<RestrictionEnzyme>, // Does not need to be saved
    backbone_lib: Vec<Backbone>,
    reading_frame: ReadingFrame,
    search_seq: Seq,
    cloning: CloningState,
}

impl Default for State {
    fn default() -> Self {
        let mut result = Self {
            ui: Default::default(),
            active: Default::default(),
            generic: vec![Default::default()],
            path_loaded: vec![Default::default()],
            portions: vec![Default::default()],
            ion_concentrations: vec![Default::default()],
            pcr: Default::default(),
            restriction_enzyme_lib: Default::default(),
            backbone_lib: Default::default(),
            reading_frame: Default::default(),
            volatile: vec![Default::default()],
            search_seq: Default::default(),
            cloning: Default::default(),
        };

        // Load the RE lib before prefs, because prefs may include loading of previously-opened files,
        // which then trigger RE match syncs.
        result.restriction_enzyme_lib = load_re_library();
        result.backbone_lib = load_backbone_library();
        result
    }
}

impl State {
    /// Add a default-settings tab, and open it.
    pub fn add_tab(&mut self) {
        self.generic.push(Default::default());
        self.ion_concentrations.push(Default::default());
        self.path_loaded.push(Default::default());
        self.portions.push(Default::default());
        self.volatile.push(Default::default());

        self.active = self.generic.len() - 1;

        // Sync items that aren't stored as part of tabs.
        self.sync_re_sites();
        self.sync_reading_frame();

        // So this tab opens on a new program run.
        self.save_prefs()
    }

    pub fn remove_tab(&mut self, i: usize) {
        let n = self.generic.len();

        if n == 1 {
            // We are closing the last tab; replace it with a blank slate.
            self.reset();
            return;
        }

        if i >= n {
            return;
        }

        self.generic.remove(i);
        self.ion_concentrations.remove(i);
        self.path_loaded.remove(i);
        self.portions.remove(i);
        self.volatile.remove(i);

        let mut tab_i_removed = None;
        for (j, tab) in self.ui.re.tabs_selected.iter().enumerate() {
            if *tab == i {
                tab_i_removed = Some(j);
            }
        }
        if let Some(j) = tab_i_removed {
            self.ui.re.tabs_selected.remove(j);
        }

        // Don't let the active tab overflow to the right; move it to the left if it would.
        // And, don't move the active tab left only if it would underflow; this effectively moves it right.
        if (self.active > 0 && self.active <= i && n > 1) || self.active + 1 >= n {
            self.active -= 1;
        }

        // So these tabs don't open on the next program run.
        self.save_prefs()
    }

    /// Convenience function, since we call this so frequently.
    pub fn get_seq(&self) -> &[Nucleotide] {
        &self.generic[self.active].seq
    }

    /// Reset data; we currently use this for making "new" data.
    pub fn reset(&mut self) {
        self.generic[self.active] = Default::default();
        self.path_loaded[self.active] = Default::default();
        self.portions[self.active] = Default::default();
        self.volatile[self.active] = Default::default();
        // todo: Ideally we reset the window title  here, but we've having trouble with variable
        // todo scope in the input function.

        self.ui.cursor_pos = None;
        self.ui.cursor_seq_i = None;
        self.ui.text_cursor_i = Some(0); // todo: For now; having trouble with cursor on empty seq
        self.ui.seq_input = String::new();
    }

    /// Load UI and related data not related to a specific sequence.
    /// This will also open all files specified in the saved UI state.
    pub fn load_prefs(&mut self, path: &Path) {
        let ui_loaded: io::Result<StateUiToSave> = load(path);

        if let Ok(ui) = ui_loaded {
            let (ui, paths_loaded) = ui.to_state();
            self.ui = ui;

            for path in &paths_loaded {
                if let Some(loaded) = load_import(path) {
                    self.load(&loaded, self.active);
                }
            }
        }
    }

    pub fn save_prefs(&self) {
        if let Err(e) = save(
            &PathBuf::from(DEFAULT_PREFS_FILE),
            &StateUiToSave::from_state(&self.ui, &self.path_loaded),
        ) {
            eprintln!("Error saving prefs: {e}");
        }
    }

    /// Runs the match search between primers and sequences. Run this when primers and sequences change.
    pub fn sync_primer_matches(&mut self, primer_i: Option<usize>) {
        let seq = &self.generic[self.active].seq.clone(); // todo; Not ideal to clone.
        let primers = match primer_i {
            Some(i) => &mut self.generic[self.active].primers[i..=i],
            // Run on all primers.
            None => &mut self.generic[self.active].primers,
        };

        for primer in primers {
            primer.volatile.matches = primer.match_to_seq(&seq);
        }
    }

    pub fn sync_pcr(&mut self) {
        self.pcr = PcrParams::new(&self.ui.pcr);
    }

    /// Identify restriction enzyme sites in the sequence.
    pub fn sync_re_sites(&mut self) {
        if self.active >= self.volatile.len() {
            eprintln!("Error: Volatile len too short for RE sync.");
            return;
        }
        self.volatile[self.active].restriction_enzyme_matches = Vec::new();

        self.volatile[self.active]
            .restriction_enzyme_matches
            .append(&mut find_re_matches(
                &self.generic[self.active].seq,
                &self.restriction_enzyme_lib,
            ));

        // This sorting aids in our up/down label alternation in the display.
        self.volatile[self.active]
            .restriction_enzyme_matches
            .sort_by(|a, b| a.seq_index.cmp(&b.seq_index));
    }

    pub fn sync_reading_frame(&mut self) {
        self.volatile[self.active].reading_frame_matches =
            find_orf_matches(self.get_seq(), self.reading_frame);
    }

    pub fn sync_search(&mut self) {
        if self.search_seq.len() >= MIN_SEARCH_LEN {
            self.volatile[self.active].search_matches =
                find_search_matches(self.get_seq(), &self.search_seq);
        } else {
            self.volatile[self.active].search_matches = Vec::new();
        }
    }

    pub fn sync_portions(&mut self) {
        for sol in &mut self.portions[self.active].solutions {
            sol.calc_amounts();
        }
    }

    pub fn sync_primer_metrics(&mut self) {
        for primer in &mut self.generic[self.active].primers {
            primer.run_calcs(&self.ion_concentrations[self.active]);
            //
            // primer.volatile[self.active].sequence_input = seq_to_str(&primer.sequence);
            //
            //
            // if let Some(metrics) = &mut primer.volatile[self.active].metrics {
            //     let dual_ended = match primer.volatile[self.active].tune_setting {
            //         TuneSetting::Both(_) => true,
            //         _ => false,
            //     };
            //
            //     metrics.update_scores(dual_ended);
            // }
            // if primer.volatile[self.active].metrics.is_none() {
            //     primer.run_calcs(&self.ion_concentrations);
            // }
        }
    }

    /// Upddate this sequence by inserting a sequence of interest. `insert_loc` uses 0-based indexing.
    pub fn insert_nucleotides(&mut self, insert: &[Nucleotide], insert_loc: usize) {
        if insert_loc == 0 {
            eprintln!("Error: Insert loc = 0");
            return;
        }

        let seq_vector = &mut self.generic[self.active].seq;

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
        for feature in &mut self.generic[self.active].features {
            if feature.range.start > insert_i {
                feature.range.start += insert.len();
            }

            if feature.range.end > insert_i {
                feature.range.end += insert.len();
            }
        }

        self.sync_seq_related(None);
    }

    /// One-based indexing. Similar to `insert_nucleotides`.
    pub fn remove_nucleotides(&mut self, range: RangeIncl) {
        let seq = &mut self.generic[self.active].seq;
        if range.end + 1 > seq.len() {
            return;
        }

        let count = range.len();

        seq.drain(range.start..=range.end);

        // Now, you have to update features affected by this insertion, shifting them left A/R.
        for feature in &mut self.generic[self.active].features {
            if feature.range.start > range.end {
                feature.range.start -= count;
            }
            if feature.range.end > range.end {
                feature.range.end -= count;
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
        self.volatile[self.active].proteins = proteins_from_seq(
            self.get_seq(),
            &self.generic[self.active].features,
            &self.volatile[self.active].cr_orf_matches,
        );

        self.ui.seq_input = seq_to_str(self.get_seq());
    }

    pub fn reset_selections(&mut self) {
        self.ui.text_selection = None;
        self.ui.selected_item = Selection::None;
    }

    /// Load state from a (our format) file.
    pub fn load(&mut self, loaded: &StateToSave, active: usize) {
        let gen = &self.generic[self.active];
        // Quick+Dirty check if we're on a new file. If so, replace it, vice adding a new tab.
        if !gen.seq.is_empty()
            || !gen.features.is_empty()
            || !gen.primers.is_empty()
            || self.path_loaded[self.active].is_some()
            || !self.portions[self.active].solutions.is_empty()
        {
            self.add_tab();
        } else {
            self.volatile.push(Default::default());
        }

        self.generic[self.active].clone_from(&loaded.generic);
        self.path_loaded[self.active].clone_from(&loaded.path_loaded);
        self.ion_concentrations[self.active].clone_from(&loaded.ion_concentrations);
        self.portions[self.active].clone_from(&loaded.portions);

        self.volatile[self.active] = Default::default();

        self.sync_pcr();
        self.sync_primer_metrics();
        self.sync_seq_related(None);
        self.ui.seq_input = seq_to_str(self.get_seq());

        self.sync_portions();
        self.reset_selections();
    }

    /// Copy the sequence of the selected text selection, feature or primer to the clipboard, if applicable.
    pub fn copy_seq(&self) {
        // Text selection takes priority.
        if let Some(selection) = &self.ui.text_selection {
            if let Some(seq) = selection.index_seq(self.get_seq()) {
                let mut ctx = ClipboardContext::new().unwrap();
                ctx.set_contents(seq_to_str(seq)).unwrap();
            }
            return;
        }

        match self.ui.selected_item {
            Selection::Feature(i) => {
                if i >= self.generic[self.active].features.len() {
                    eprintln!("Invalid feature in selection");
                    return;
                }
                let feature = &self.generic[self.active].features[i];
                if let Some(seq) = feature.range.index_seq(self.get_seq()) {
                    let mut ctx = ClipboardContext::new().unwrap();
                    ctx.set_contents(seq_to_str(seq)).unwrap();
                }
            }
            Selection::Primer(i) => {
                let primer = &self.generic[self.active].primers[i];

                let mut ctx = ClipboardContext::new().unwrap();
                ctx.set_contents(seq_to_str(&primer.sequence)).unwrap();
            }
            _ => (),
        }
    }
}

fn main() {
    let mut state = State::default();

    state.load_prefs(&PathBuf::from_str(DEFAULT_PREFS_FILE).unwrap());

    // Initial load  hierarchy:
    // - Path argument (e.g. file association)
    // - Last opened files
    // - Quicksave

    let mut loaded_from_arg = false;
    let (path, window_title_initial) = {
        let mut p = PathBuf::from_str(QUICKSAVE_FILE).unwrap();

        // Windows and possibly other operating systems, if attempting to use your program to natively
        // open a file type, will use command line arguments to indicate this. Determine if the program
        // is being launched this way, and if so, open the file.
        let args: Vec<String> = env::args().collect();
        if args.len() > 1 {
            let temp = &args[1];
            p = PathBuf::from_str(temp).unwrap();
            loaded_from_arg = true;
        }

        let window_title = get_window_title(&p);
        (p, window_title)
    };

    let mut prev_paths_loaded = false;
    for path in &state.path_loaded {
        if path.is_some() {
            prev_paths_loaded = true;
        }
    }

    // Load from the argument or quicksave A/R.
    if loaded_from_arg || !prev_paths_loaded {
        if let Some(loaded) = load_import(&path) {
            println!("LOADING QS");
            state.load(&loaded, state.active);
        }
    }

    state.sync_seq_related(None);

    state.reset_selections();

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
