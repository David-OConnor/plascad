mod primer;

use eframe::{
    egui,
    egui::{Color32, Context, Ui},
};

use crate::{
    primer::{Primer, PrimerMetrics},
    seq_from_str, State,
};

pub const WINDOW_WIDTH: f32 = 1200.;
pub const WINDOW_HEIGHT: f32 = 800.;

pub const WINDOW_TITLE: &str = "Plasmid check";

pub const ROW_SPACING: f32 = 22.;
pub const COL_SPACING: f32 = 30.;

#[derive(Clone, Copy)]
pub enum Page {
    Primers,
    // Sequence,
    // Enzymes,
    // Features,
}

impl Default for Page {
    fn default() -> Self {
        Self::Primers
    }
}

// // todo: Move A/R
// /// Used to determine which side of a primer we can extend or remove from in order to optimize it.
// #[derive(Clone, Copy, PartialEq)]
// enum TunableEnd {
//     None,
//     Left,
//     Right,
//     Both,
// }

// impl Default for TunableEnd {
//     fn default() -> Self {
//         Self::None
//     }
// }

#[derive(Default)]
pub struct PrimerTableCol {
    primer: Primer,
    /// Editing is handled using this string; we convert the string to our nucleotide sequence as needed.
    sequence_input: String,
    description: String,
    metrics: Option<PrimerMetrics>,
    // tunable_end: TunableEnd,
    /// These fields control if a given primer end is fixed (Eg marking the start of an insert,
    /// marking the insert point in a vector etc) or if we can tune its length to optimize the primer.
    tunable_5p: bool,
    tunable_3p: bool,
}

impl PrimerTableCol {
    /// Perform calculations on primer quality and related data. Run this when the sequence changes.
    pub fn run_calcs(&mut self) {
        self.primer.sequence = seq_from_str(&self.sequence_input);
        self.metrics = self.primer.calc_metrics();
    }
}

fn page_selector(state: &mut State, ui: &mut Ui) {
    ui.add_space(ROW_SPACING);
}

pub fn draw(state: &mut State, ctx: &Context) {
    egui::CentralPanel::default().show(ctx, |ui| {
        let mut visuals = ctx.style().visuals.clone();
        // visuals.override_text_color = Some(Color32::from_rgb(255, 0, 0));
        visuals.override_text_color = Some(Color32::LIGHT_GRAY);
        ctx.set_visuals(visuals);

        egui::containers::ScrollArea::vertical().show(ui, |ui| match state.ui.page {
            Page::Primers => primer::primer_page(state, ui),
        });
    });
}
