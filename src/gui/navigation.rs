//! This module contains code related to navigation buttons.

use std::fmt::Display;

use eframe::egui::{Color32, RichText, Ui};

use crate::{gui::COL_SPACING, State};

pub const NAV_BUTTON_COLOR: Color32 = Color32::from_rgb(0, 00, 110);

#[derive(Clone, Copy, PartialEq)]
pub enum Page {
    /// Primer design and QC, including for cloning
    /// (Replacement name: Sequence?
    Sequence,
    /// A circular "graphical map" of the plasmid
    Map,
    Features,
    Primers,
    /// Determine optimal PCR parameters
    Pcr,
    Portions,
    // Enzymes,
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
            Self::Pcr => "PCR",
            Self::Features => "Features",
            Self::Primers => "Primers",
            Self::Portions => "Mixing portions",
        }
        .to_owned();
        write!(f, "{}", str)
    }
}

pub fn page_selector(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        page_button(&mut state.ui.page, Page::Sequence, ui);
        page_button(&mut state.ui.page, Page::Map, ui);
        page_button(&mut state.ui.page, Page::Features, ui);
        page_button(&mut state.ui.page, Page::Primers, ui);
        page_button(&mut state.ui.page, Page::Pcr, ui);
        // page_button(&mut state.ui.page, Page::Portions, ui);
    });
}

fn page_button<T: PartialEq + ToString>(page_state: &mut T, page: T, ui: &mut Ui) {
    let color = if *page_state == page {
        Color32::GREEN
    } else {
        Color32::WHITE
    };

    if ui
        .button(
            RichText::new(page.to_string())
                .color(color)
                .background_color(NAV_BUTTON_COLOR),
        )
        .clicked()
    {
        *page_state = page;
    }

    ui.add_space(COL_SPACING / 2.);
}

/// This is used for selecting what is displayed in the sequence view, ie view or edit.
#[derive(Clone, Copy, PartialEq)]
pub enum PageSeq {
    EditSeq,
    EditSlic,
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
            Self::EditSeq => "Edit sequence",
            Self::EditSlic => "SLIC/FC cloning",
            Self::View => "View sequence",
        }
        .to_owned();
        write!(f, "{}", str)
    }
}

pub fn page_seq_selector(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        page_button(&mut state.ui.page_seq, PageSeq::EditSeq, ui);
        page_button(&mut state.ui.page_seq, PageSeq::View, ui);
        page_button(&mut state.ui.page_seq, PageSeq::EditSlic, ui);
    });
}

/// This is used for selecting what is displayed above the sequence view, ie various tabular editors.
#[derive(Clone, Copy, PartialEq)]
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
        page_button(&mut state.ui.page_seq_top, PageSeqTop::Primers, ui);
        page_button(&mut state.ui.page_seq_top, PageSeqTop::Features, ui);
        page_button(&mut state.ui.page_seq_top, PageSeqTop::None, ui);
    });
}
