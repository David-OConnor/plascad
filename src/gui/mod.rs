use std::str::FromStr;

use eframe::{
    egui,
    egui::{Color32, Context, Key, ScrollArea, TextEdit, Ui},
};
use navigation::Page;

use crate::{
    file_io::save::{save, StateToSave, DEFAULT_SAVE_FILE},
    gui::primer_qc::primer_details,
    State,
};

mod circle;
mod feature_overlay;
mod features;
mod metadata;
pub mod navigation;
mod pcr;
mod portions;
mod primer_arrow;
pub mod primer_qc;
mod save;
pub mod seq_view;
pub mod sequence;
// pub for a few consts

pub const WINDOW_WIDTH: f32 = 1300.;
pub const WINDOW_HEIGHT: f32 = 1_000.;

pub const WINDOW_TITLE: &str = "PlasCAD";

pub const ROW_SPACING: f32 = 22.;
pub const COL_SPACING: f32 = 30.;

pub fn int_field(val: &mut usize, label: &str, ui: &mut Ui) {
    ui.label(label);
    let mut entry = val.to_string();
    let response = ui.add(TextEdit::singleline(&mut entry).desired_width(40.));
    if response.changed() {
        *val = entry.parse().unwrap_or(0);
    }
}

/// Get a text-representation of the cursor index; a slightly processed version of the raw index.
/// We use this on the sequence and circle views.
pub fn get_cursor_text(cursor_seq_i: Option<usize>, seq_len: usize) -> String {
    match cursor_seq_i {
        Some(p) => {
            if p + 1 <= seq_len {
                // + 1, as the convention is to use 1-based indexing vice 0.
                (p + 1).to_string()
                // This occurs if the cursor is on the last row, right of the last NT.
            } else {
                String::new()
            }
        }
        None => String::new(),
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

            ui.label("Name: ");
            ui.add(TextEdit::singleline(&mut state.plasmid_name).desired_width(140.));

            ui.add_space(COL_SPACING);

            save::save_section(state, ui);
        });

        ui.add_space(ROW_SPACING);

        ScrollArea::vertical().show(ui, |ui| match state.ui.page {
            Page::Sequence => sequence::seq_page(state, ui),
            Page::Map => circle::circle_page(state, ui),
            Page::Features => features::features_page(state, ui),
            Page::Primers => primer_details(state, ui),
            Page::Pcr => pcr::pcr_page(state, ui),
            Page::Metadata => metadata::metadata_page(state, ui),
            _ => (),
            // Page::Portions => portions::portions_page(state, ui),
        });
    });
}
