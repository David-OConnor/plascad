//! Data structures relalted to general state.

use std::{
    io,
    path::{Path, PathBuf},
    time::Instant,
};

use copypasta::{ClipboardContext, ClipboardProvider};
use eframe::egui::Context;
use na_seq::{
    insert_into_seq,
    ligation::LigationFragment,
    re_lib::load_re_library,
    restriction_enzyme::{find_re_matches, ReMatch, RestrictionEnzyme},
    seq_to_str_lower, Nucleotide, Seq,
};

use crate::{
    ab1::SeqRecordAb1,
    alignment::AlignmentState,
    backbones::{load_backbone_library, Backbone},
    cloning::CloningState,
    file_io::{
        save::{load, load_import, save, PrefsToSave, StateToSave, DEFAULT_PREFS_FILE},
        GenericData,
    },
    gui,
    gui::navigation::Tab,
    misc_types::{find_search_matches, SearchMatch, MIN_SEARCH_LEN},
    pcr::PcrParams,
    portions::PortionsState,
    primer::IonConcentrations,
    protein::{proteins_from_seq, sync_cr_orf_matches, Protein},
    reading_frame::{find_orf_matches, ReadingFrame, ReadingFrameMatch},
    tags::TagMatch,
    util::RangeIncl,
    Selection, StateUi, PREFS_SAVE_INTERVAL,
};

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

/// Note: use of serde traits here and on various sub-structs are for saving and loading.
pub struct State {
    pub ui: StateUi,
    /// Ie tab, file etc.
    pub active: usize,
    // todo: Consider grouping generic, path_loaded, portions, and similar in a single vec.
    // todo: Do that after your initial tab approach works.
    /// Data that is the most fundamental to persistent state, and shared between save formats.
    /// Index corresponds to `active`.
    pub generic: Vec<GenericData>,
    /// Index corresponds to `active`.
    pub ab1_data: Vec<SeqRecordAb1>,
    // Used to determine which file to save to, if applicable.
    // file_active: Option<Tab>,
    /// Index corresponds to `active`.
    pub tabs_open: Vec<Tab>,
    /// Index corresponds to `active`.
    pub portions: Vec<PortionsState>,
    /// Index corresponds to `active`.
    pub volatile: Vec<StateVolatile>,
    /// Used for PCR.
    // todo: YOu may need to go back to per-tab ion concentrations.
    // ion_concentrations: Vec<IonConcentrations>,
    pub ion_concentrations: IonConcentrations,
    pub pcr: PcrParams,
    pub restriction_enzyme_lib: Vec<RestrictionEnzyme>, // Does not need to be saved
    pub backbone_lib: Vec<Backbone>,
    pub reading_frame: ReadingFrame,
    pub search_seq: Seq,
    pub cloning: CloningState,
    pub alignment: AlignmentState,
}

impl Default for State {
    fn default() -> Self {
        let mut result = Self {
            ui: Default::default(),
            active: Default::default(),
            generic: vec![Default::default()],
            ab1_data: vec![Default::default()],
            tabs_open: vec![Default::default()],
            portions: vec![Default::default()],
            // ion_concentrations: vec![Default::default()],
            ion_concentrations: Default::default(),
            pcr: Default::default(),
            restriction_enzyme_lib: Default::default(),
            backbone_lib: Default::default(),
            reading_frame: Default::default(),
            volatile: vec![Default::default()],
            search_seq: Default::default(),
            cloning: Default::default(),
            alignment: Default::default(),
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
        // self.ion_concentrations.push(Default::default());
        // self.tabs_open.push(Default::default());
        self.portions.push(Default::default());
        self.volatile.push(Default::default());
        self.ab1_data.push(Default::default());

        self.active = self.generic.len() - 1;

        // todo: DRY with reset
        self.ui.cursor_pos = None;
        self.ui.cursor_seq_i = None;
        self.ui.text_cursor_i = Some(0);
        self.ui.seq_input = String::new();

        // Sync items that aren't stored as part of tabs.
        self.sync_re_sites();
        self.sync_reading_frame();

        // println!("Saving prefs");
        // self.save_prefs()
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
        self.ab1_data.remove(i);
        // self.ion_concentrations.remove(i);
        self.tabs_open.remove(i);
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
        self.tabs_open[self.active] = Default::default();
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
    /// This will also open all files specified in the saved preferences.
    pub fn load_prefs(&mut self, path: &Path) {
        let prefs_loaded: io::Result<PrefsToSave> = load(path);

        if let Ok(prefs) = prefs_loaded {
            let (ui, tabs_open, ion_concentrations) = prefs.to_state();
            self.ui = ui;
            self.ion_concentrations = ion_concentrations;

            for tab in &tabs_open {
                if let Some(path) = &tab.path {
                    if let Some(loaded) = load_import(path) {
                        self.load(&loaded);
                    }
                }
            }
        }
    }

    pub fn save_prefs(&self) {
        if let Err(e) = save(
            &PathBuf::from(DEFAULT_PREFS_FILE),
            &PrefsToSave::from_state(&self.ui, &self.tabs_open, &self.ion_concentrations),
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
            // primer.run_calcs(&self.ion_concentration[self.active]);
            primer.run_calcs(&self.ion_concentrations);
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

    /// Upddate this sequence by inserting a sequence of interest. Shifts features based on the insert.
    pub fn insert_nucleotides(&mut self, insert: &[Nucleotide], insert_loc: usize) {
        insert_into_seq(&mut self.generic[self.active].seq, insert, insert_loc).ok();

        let insert_i = insert_loc - 1; // 1-based indexing.

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

        // We have removed RE sync here for now, because it is slow. Call  when able.
        // self.sync_re_sites();

        self.sync_reading_frame();
        self.sync_search();

        sync_cr_orf_matches(self);

        self.volatile[self.active].proteins = proteins_from_seq(
            self.get_seq(),
            &self.generic[self.active].features,
            &self.volatile[self.active].cr_orf_matches,
        );

        self.ui.seq_input = seq_to_str_lower(self.get_seq());
    }

    pub fn reset_selections(&mut self) {
        self.ui.text_selection = None;
        self.ui.selected_item = Selection::None;
    }

    /// Load a single tab from file.
    pub fn load(&mut self, loaded: &StateToSave) {
        let gen = &self.generic[self.active];

        // Check if the current tab is a new file. If so, don't add a new tab.
        if !gen.seq.is_empty()
            || !gen.features.is_empty()
            || !gen.primers.is_empty()
            || !self.portions[self.active].solutions.is_empty()
            || !self.ab1_data[self.active].sequence.is_empty()
        {
            self.add_tab();

            self.tabs_open.push(Default::default());
        }

        // todo: Fix this for AB1 loading.
        let ab1 = !loaded.ab1_data.sequence.is_empty();
        let tab = Tab {
            path: loaded.path_loaded.clone(),
            ab1,
        };

        self.generic[self.active].clone_from(&loaded.generic);
        // self.ion_concentrations[self.active].clone_from(&loaded.ion_concentrations);
        self.portions[self.active].clone_from(&loaded.portions);
        self.ab1_data[self.active].clone_from(&loaded.ab1_data);
        self.tabs_open[self.active] = tab;

        self.volatile[self.active] = Default::default();

        self.sync_pcr();
        self.sync_primer_metrics();
        self.sync_seq_related(None);
        self.sync_re_sites();

        self.ui.seq_input = seq_to_str_lower(self.get_seq());

        self.sync_portions();
        self.reset_selections();

        // So these tabs opens on a new program run.
        self.save_prefs(); // Save opened tabs.
    }

    /// Copy the sequence of the selected text selection, feature or primer to the clipboard, if applicable.
    pub fn copy_seq(&self) {
        // Text selection takes priority.
        if let Some(selection) = &self.ui.text_selection {
            if let Some(seq) = selection.index_seq(self.get_seq()) {
                let mut ctx = ClipboardContext::new().unwrap();
                ctx.set_contents(seq_to_str_lower(seq)).unwrap();
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
                    ctx.set_contents(seq_to_str_lower(seq)).unwrap();
                }
            }
            Selection::Primer(i) => {
                let primer = &self.generic[self.active].primers[i];

                let mut ctx = ClipboardContext::new().unwrap();
                ctx.set_contents(seq_to_str_lower(&primer.sequence))
                    .unwrap();
            }
            _ => (),
        }
    }
}

/// This struct contains state that does not need to persist between sessesions or saves, but is not
/// a good fit for `StateUi`. This is, generally, calculated data from persistent staet.
#[derive(Default)]
pub struct StateVolatile {
    pub restriction_enzyme_matches: Vec<ReMatch>,
    pub re_digestion_products: Vec<LigationFragment>,
    pub reading_frame_matches: Vec<ReadingFrameMatch>,
    pub tag_matches: Vec<TagMatch>,
    pub search_matches: Vec<SearchMatch>,
    /// Used for automatically determining which reading frame to use, and the full frame,
    /// for a given coding-region feature.
    pub cr_orf_matches: Vec<(usize, ReadingFrameMatch)>,
    pub proteins: Vec<Protein>,
}
