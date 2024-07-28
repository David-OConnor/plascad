mod pcr;
mod portions;
pub mod primer;
pub mod seq_view; // pub for a few consts

use std::{fs::OpenOptions, io};

use bio::io::fasta;
use eframe::{
    egui,
    egui::{Color32, Context, Key, RichText, ScrollArea, Ui},
};
use egui_file::FileDialog;

use crate::{
    save::{export_fasta, import_fasta, save, StateToSave, DEFAULT_FASTA_FILE, DEFAULT_SAVE_FILE},
    sequence::make_seq_str,
    State,
};

pub const WINDOW_WIDTH: f32 = 1300.;
pub const WINDOW_HEIGHT: f32 = 800.;

pub const WINDOW_TITLE: &str = "Plasmid tools";

pub const ROW_SPACING: f32 = 22.;
pub const COL_SPACING: f32 = 30.;

#[derive(Clone, Copy, PartialEq)]
pub enum Page {
    /// Primer design and QC, including for cloning
    Primers,
    /// Determine optimal PCR parameters
    Pcr,
    Portions,
    // Sequence,
    // Enzymes,
    // Features,
}

impl Default for Page {
    fn default() -> Self {
        Self::Primers
    }
}

impl Page {
    pub fn to_str(self) -> String {
        match self {
            Self::Primers => "Primers",
            Self::Pcr => "PCR",
            Self::Portions => "Mixing portions",
        }
        .to_owned()
    }
}

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

impl PageSeq {
    pub fn to_str(self) -> String {
        match self {
            Self::EditSeq => "Edit squence",
            Self::EditSlic => "Edit SLIC/FC squences",
            Self::View => "View sequence",
        }
        .to_owned()
    }
}

fn page_button(page_state: &mut Page, page: Page, ui: &mut Ui) {
    let color = if *page_state == page {
        Color32::GREEN
    } else {
        Color32::WHITE
    };

    if ui
        .button(RichText::new(page.to_str()).color(color))
        .clicked()
    {
        *page_state = page;
    }

    ui.add_space(COL_SPACING);
}

fn page_selector(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        page_button(&mut state.ui.page, Page::Primers, ui);
        page_button(&mut state.ui.page, Page::Pcr, ui);
        page_button(&mut state.ui.page, Page::Portions, ui);
    });
}

// todo: Use to_string and partialEq traits instead of duplicating the other page.
fn page_seq_button(page_state: &mut PageSeq, page: PageSeq, ui: &mut Ui) {
    let color = if *page_state == page {
        Color32::GREEN
    } else {
        Color32::WHITE
    };

    if ui
        .button(RichText::new(page.to_str()).color(color))
        .clicked()
    {
        *page_state = page;
    }

    ui.add_space(COL_SPACING / 2.);
}

pub fn page_seq_selector(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        page_seq_button(&mut state.ui.page_seq, PageSeq::EditSeq, ui);
        page_seq_button(&mut state.ui.page_seq, PageSeq::EditSlic, ui);
        page_seq_button(&mut state.ui.page_seq, PageSeq::View, ui);
    });
}

/// Ui elements for saving and loading data in various file formats. This includes our own format,
/// FASTA, and (eventually) SnapGene's DNA format.
fn save_section(state: &mut State, ui: &mut Ui) {
    if ui
        .button("Save")
        .on_hover_text("Save primer data. Ctrl + S")
        .clicked()
    {
        if let Err(e) = save(DEFAULT_SAVE_FILE, &StateToSave::from_state(state)) {
            println!("Error saving: {e}");
        }
    }

    if ui
        .button("Import FASTA")
        .on_hover_text("Import a sequence in the FASTA format")
        .clicked()
    {
        let mut dialog = FileDialog::open_file(state.ui.opened_file.clone());
        dialog.open();
        state.ui.open_file_dialog_import = Some(dialog);
    }

    if ui
        .button("Export FASTA")
        .on_hover_text("Export the sequence in the FASTA format")
        .clicked()
    {
        let mut dialog = FileDialog::save_file(state.ui.opened_file.clone())
            .default_filename(DEFAULT_FASTA_FILE);
        dialog.open();
        state.ui.open_file_dialog_export = Some(dialog);
    }

    if let Some(dialog) = &mut state.ui.open_file_dialog_import {
        if dialog.show(ui.ctx()).selected() {
            if let Some(path) = dialog.path() {
                state.ui.opened_file = Some(path.to_owned());

                if let Ok(seq) = import_fasta(path) {
                    state.seq = seq;
                    state.ui.seq_input = make_seq_str(&state.seq);
                    state.sync_re_sites();
                }

                state.ui.opened_file = None;
                state.ui.open_file_dialog_import = None;
            }
        }
    } else if let Some(dialog) = &mut state.ui.open_file_dialog_export {
        if dialog.show(ui.ctx()).selected() {
            if let Some(path) = dialog.path() {
                state.ui.opened_file = Some(path.to_owned());

                if let Err(e) = export_fasta(&state.seq, path) {
                    eprintln!("Error exporting to FASTA: {:?}", e);
                };

                state.ui.opened_file = None;
                state.ui.open_file_dialog_export = None;
            }
        }
    }
}

pub fn draw(state: &mut State, ctx: &Context) {
    ctx.input(|ip| {
        if ip.key_pressed(Key::A) && ip.modifiers.ctrl {
            state.primer_data.push(Default::default());
        }

        if ip.key_pressed(Key::S) && ip.modifiers.ctrl {
            if let Err(e) = save(DEFAULT_SAVE_FILE, &StateToSave::from_state(state)) {
                println!("Error saving: {e}");
            }
        }

        state.ui.cursor_pos = ip.pointer.hover_pos().map(|pos| (pos.x, pos.y));
    });

    egui::CentralPanel::default().show(ctx, |ui| {
        // todo: This section DRY with seq viewx.

        let mut visuals = ctx.style().visuals.clone();
        // visuals.override_text_color = Some(Color32::from_rgb(255, 0, 0));
        visuals.override_text_color = Some(Color32::LIGHT_GRAY);
        ctx.set_visuals(visuals);

        ui.horizontal(|ui| {
            page_selector(state, ui);

            ui.add_space(COL_SPACING);

            save_section(state, ui);
        });

        ui.add_space(ROW_SPACING);

        ScrollArea::vertical().show(ui, |ui| match state.ui.page {
            Page::Primers => primer::primer_page(state, ui),
            Page::Pcr => pcr::pcr_page(state, ui),
            Page::Portions => portions::portions_page(state, ui),
        });
    });
}
