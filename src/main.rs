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
use std::{env, path::PathBuf, str::FromStr, sync::Arc};

use bincode::{Decode, Encode};
use bio::alignment::Alignment;
use cloning::{CloneStatus, CloningInsertData};
use copypasta::ClipboardProvider;
use eframe::{
    self,
    egui::{self},
};
use egui_file_dialog::{FileDialog, FileDialogConfig};
use file_io::save::{load_import, QUICKSAVE_FILE};
use gui::navigation::{Page, PageSeq};
use na_seq::{
    ligation::LigationFragment,
    restriction_enzyme::{ReMatch, RestrictionEnzyme},
    AaIdent, AminoAcid, Nucleotide, Seq,
};
use protein::Protein;
use reading_frame::ReadingFrameMatch;
use state::State;

use crate::{
    backbones::{Backbone, BackboneFilters},
    cloning::BackboneSelected,
    file_io::{
        save::{DEFAULT_DNA_FILE, DEFAULT_FASTA_FILE, DEFAULT_GENBANK_FILE, DEFAULT_PREFS_FILE},
        FileDialogs, GenericData,
    },
    gui::{navigation::PageSeqTop, WINDOW_HEIGHT, WINDOW_WIDTH},
    misc_types::{FeatureDirection, FeatureType, SearchMatch},
    pcr::{PcrUi, PolymeraseType},
    primer::{Primer, TM_TARGET},
    tags::TagMatch,
    util::{get_window_title, RangeIncl},
};

mod ab1;
mod alignment;
mod alignment_map;
mod backbones;
mod cloning;
mod external_websites;
mod feature_db_load;
mod file_io;
mod gui;
mod melting_temp_calcs;
mod misc_types;
mod pcr;
mod portions;
mod primer;
mod primer_metrics;
mod protein;
mod reading_frame;
mod save_compat;
mod solution_helper;
mod state;
mod tags;
mod toxic_proteins;
mod util;

type Color = (u8, u8, u8); // RGB

// todo: Eventually, implement a system that automatically checks for changes, and don't
// todo save to disk if there are no changes.
const PREFS_SAVE_INTERVAL: u64 = 60; // Save user preferences this often, in seconds.

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

/// UI state for restriction enzymes.
pub struct ReUi {
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
    seq_edit_lock: bool,
    ab1_start_i: usize,
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
            aa_ident_disp: AaIdent::OneLetter,
            pdb_error_received: false,
            re: Default::default(),
            backbone_filters: Default::default(),
            seq_edit_lock: true,
            ab1_start_i: Default::default(),
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

fn main() {
    let mut state = State::default();

    // todo: Temp to test BAM
    let am = alignment_map::import(&PathBuf::from_str("../../Desktop/test.bam").unwrap()).unwrap();
    println!("Alignment map loaded: {:?}", am);

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
    for tab in &state.tabs_open {
        if tab.path.is_some() {
            prev_paths_loaded = true;
        }
    }

    // todo: Consider a standalone method for loaded-from-arg,

    // Load from the argument or quicksave A/R.
    if loaded_from_arg || !prev_paths_loaded {
        println!("Loading from quicksave or arg: {:?}", path);
        if let Some(loaded) = load_import(&path) {
            state.load(&loaded);
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
        ..Default::default()
    };

    eframe::run_native(
        &window_title_initial,
        options,
        Box::new(|_cc| Ok(Box::new(state))),
    )
    .unwrap();
}
