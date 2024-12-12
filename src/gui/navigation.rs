//! This module contains code related to navigation buttons.

use std::{fmt::Display, path::PathBuf};

use bincode::{Decode, Encode};
use eframe::egui::{Color32, RichText, Ui};
use na_seq::seq_to_str_lower;

use crate::{
    gui::{select_color_text, set_window_title, COL_SPACING, ROW_SPACING},
    state::State,
    util::name_from_path,
};

pub const NAV_BUTTON_COLOR: Color32 = Color32::from_rgb(0, 0, 110);
pub const TAB_BUTTON_COLOR: Color32 = Color32::from_rgb(40, 80, 110);
pub const DEFAULT_TAB_NAME: &str = "New plasmid";

// todo: Alt name: Path loaded
#[derive(Encode, Decode, Clone, Default, Debug)]
pub struct Tab {
    pub path: Option<PathBuf>,
    pub ab1: bool, // todo: Enum if you add a third category.
}

/// Used in several GUI components to get data from open tabs.
/// Note: For name, we currently default to file name (with extension), then
/// plasmid name, then a default. See if you want to default to plasmid name.
///
/// Returns the name, and the tab index.
pub fn get_tab_names(
    tabs: &[Tab],
    plasmid_names: &[&str],
    abbrev_name: bool,
) -> Vec<(String, usize)> {
    let mut result = Vec::new();

    for (i, p) in tabs.iter().enumerate() {
        let name = name_from_path(&p.path, plasmid_names[i], abbrev_name);
        result.push((name, i));
    }

    result
}

#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum Page {
    /// Primer design and QC, including for cloning
    Sequence,
    /// A circular "graphical map" of the plasmid
    Map,
    Features,
    Primers,
    Proteins,
    /// Determine optimal PCR parameters
    Pcr,
    Alignment,
    Portions,
    Metadata,
    Ligation,
    /// A simplified cloning process, for both PCR and RE-based cloning.
    Cloning,
    /// i.e. Sanger sequencing data. This page is fundamentally different from the others;
    /// it is only for .ab1 files, and is selected automatically, vice from the menu.
    Ab1,
}

impl Default for Page {
    fn default() -> Self {
        Self::Sequence
    }
}

impl Display for Page {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            Self::Sequence => "Sequence",
            Self::Map => "Map",
            Self::Features => "Features",
            Self::Primers => "Primers",
            Self::Proteins => "Proteins",
            Self::Pcr => "PCR",
            Self::Alignment => "Align",
            Self::Portions => "Mixing",
            Self::Metadata => "Data",
            Self::Ligation => "Digest",
            Self::Cloning => "Clone",
            Self::Ab1 => "AB1",
        }
        .to_owned();
        write!(f, "{}", str)
    }
}

/// Selects which tab (ie file) is active
pub fn tab_selector(state: &mut State, ui: &mut Ui) {
    // Note: This assumes generic, paths_loaded etc are always the same length. (Which they *should* be.)

    let mut tab_removed = None;

    ui.horizontal(|ui| {
        if ui
            .button("New")
            .on_hover_text("Create and open a new file. (Ctrl + N)")
            .clicked()
        {
            state.add_tab();
            state.tabs_open.push(Default::default());
        }
        ui.add_space(COL_SPACING);

        let plasmid_names: &Vec<_> = &state
            .generic
            .iter()
            .map(|v| v.metadata.plasmid_name.as_str())
            .collect();

        for (name, i) in get_tab_names(&state.tabs_open, plasmid_names, false) {
            // todo: DRY with page selectors.

            let button = ui.button(
                select_color_text(&name, i == state.active).background_color(TAB_BUTTON_COLOR),
            );

            if button.clicked() {
                state.active = i;
                set_window_title(&state.tabs_open[i], ui);

                // todo: Apt state sync fn for this?
                state.ui.seq_input = seq_to_str_lower(state.get_seq()); // todo: Move seq_input to an indexed vector?
                                                                        // todo: Cache these instead?
                                                                        // state.sync_seq_related(None);
            }

            if button.middle_clicked() {
                tab_removed = Some(i);
            }

            ui.add_space(COL_SPACING / 2.);
        }

        // todo: Right-align?
        ui.add_space(2. * ROW_SPACING);

        if ui
            .button(
                RichText::new("Close active tab")
                    .color(Color32::WHITE)
                    .background_color(Color32::DARK_RED),
            )
            .on_hover_text("Shortcut: Middle click the tab to close it.")
            .clicked
        {
            tab_removed = Some(state.active)
        };
    });

    if let Some(i) = tab_removed {
        state.remove_tab(i);
    }
}

pub fn page_selector(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        page_button(&mut state.ui.page, Page::Sequence, ui, true);
        page_button(&mut state.ui.page, Page::Map, ui, true);
        page_button(&mut state.ui.page, Page::Features, ui, true);
        page_button(&mut state.ui.page, Page::Primers, ui, true);
        page_button(&mut state.ui.page, Page::Proteins, ui, true);
        page_button(&mut state.ui.page, Page::Cloning, ui, true);
        page_button(&mut state.ui.page, Page::Pcr, ui, true);
        page_button(&mut state.ui.page, Page::Alignment, ui, true);
        page_button(&mut state.ui.page, Page::Ligation, ui, true);
        page_button(&mut state.ui.page, Page::Metadata, ui, true);
        page_button(&mut state.ui.page, Page::Portions, ui, true);
    });
}

pub fn page_button<T: PartialEq + ToString>(page_state: &mut T, page: T, ui: &mut Ui, space: bool) {
    if ui
        .button(
            select_color_text(&page.to_string(), *page_state == page)
                .background_color(NAV_BUTTON_COLOR),
        )
        .clicked()
    {
        *page_state = page;
    }

    if space {
        ui.add_space(COL_SPACING / 2.);
    }
}

/// This is used for selecting what is displayed in the sequence view, ie view or edit.
#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum PageSeq {
    EditRaw,
    // EditSlic,
    View,
}

impl Default for PageSeq {
    fn default() -> Self {
        Self::View
    }
}

impl Display for PageSeq {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            Self::EditRaw => "Edit raw",
            // Self::EditSlic => "SLIC/FC cloning",
            Self::View => "Sequence",
        }
        .to_owned();
        write!(f, "{}", str)
    }
}

pub fn page_seq_selector(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        page_button(&mut state.ui.page_seq, PageSeq::EditRaw, ui, true);
        page_button(&mut state.ui.page_seq, PageSeq::View, ui, true);
    });
}

/// This is used for selecting what is displayed above the sequence view, ie various tabular editors.
#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum PageSeqTop {
    Primers,
    Features,
    None,
}

impl Default for PageSeqTop {
    fn default() -> Self {
        Self::None
    }
}

impl Display for PageSeqTop {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            Self::Primers => "Primers",
            Self::Features => "Features",
            Self::None => "None",
        }
        .to_owned();
        write!(f, "{}", str)
    }
}

pub fn page_seq_top_selector(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        ui.label("Display above sequence:");

        page_button(&mut state.ui.page_seq_top, PageSeqTop::None, ui, true);
        page_button(&mut state.ui.page_seq_top, PageSeqTop::Features, ui, true);
        page_button(&mut state.ui.page_seq_top, PageSeqTop::Primers, ui, true);
    });
}
