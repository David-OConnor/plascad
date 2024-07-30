mod features;
pub mod navigation;
mod pcr;
mod portions;
pub mod primer;
pub mod seq_view;
// pub for a few consts

use std::fmt::Display;

use eframe::{
    egui,
    egui::{Color32, Context, Key, ScrollArea, TextEdit, Ui},
};
use egui_file::FileDialog;
use navigation::Page;

use crate::{
    save::{export_fasta, import_fasta, save, StateToSave, DEFAULT_FASTA_FILE, DEFAULT_SAVE_FILE},
    sequence::seq_to_str,
    State,
};

pub const WINDOW_WIDTH: f32 = 1300.;
pub const WINDOW_HEIGHT: f32 = 800.;

pub const WINDOW_TITLE: &str = "Plasmid tools";

pub const ROW_SPACING: f32 = 22.;
pub const COL_SPACING: f32 = 30.;

pub fn int_field(val: &mut usize, label: &str, ui: &mut Ui) {
    ui.label("Start:");
    let mut entry = val.to_string();
    let response = ui.add(TextEdit::singleline(&mut entry).desired_width(40.));
    if response.changed() {
        *val = entry.parse().unwrap_or(0);
    }
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
                    state.ui.seq_input = seq_to_str(&state.seq);
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
            navigation::page_selector(state, ui);

            ui.add_space(COL_SPACING);

            save_section(state, ui);
        });

        ui.add_space(ROW_SPACING);

        ScrollArea::vertical().show(ui, |ui| match state.ui.page {
            Page::Sequence => primer::primer_page(state, ui),
            Page::Pcr => pcr::pcr_page(state, ui),
            Page::Portions => portions::portions_page(state, ui),
        });
    });
}
