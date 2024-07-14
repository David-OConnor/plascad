pub mod primer;

use eframe::{
    egui,
    egui::{Color32, Context, Key, Ui},
};

use crate::State;

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

fn page_selector(state: &mut State, ui: &mut Ui) {
    ui.add_space(ROW_SPACING);
}

pub fn draw(state: &mut State, ctx: &Context) {
    let input = ctx.input(|ip| {
        if ip.key_pressed(Key::A) && ip.modifiers.ctrl {
            state.ui.primer_cols.push(Default::default());
        }

        if ip.key_pressed(Key::S) && ip.modifiers.ctrl {
            // save()
        }
    });

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
