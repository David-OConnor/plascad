//! UI page for mixing portions. (growth media, stock solutions etc)

use eframe::egui::Ui;

use crate::{gui::ROW_SPACING, State};

pub fn portions_page(state: &mut State, ui: &mut Ui) {
    ui.heading("Mixing portions");

    ui.add_space(ROW_SPACING);

    ui.heading("Growth media");
}
