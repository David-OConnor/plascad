//! Contains code for viewing AB1 sequencing data, e.g. from Sanger sequencing.

use eframe::egui::Ui;

use crate::ab1::SeqRecordAb1;

pub fn ab1_page(data: &SeqRecordAb1, ui: &mut Ui) {
    ui.heading("AB1 sequencing view");
}
