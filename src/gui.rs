use eframe::{
    egui,
    egui::{Color32, Context, RichText},
};

use crate::State;

pub const WINDOW_WIDTH: f32 = 1100.;
pub const WINDOW_HEIGHT: f32 = 1_000.;

pub const WINDOW_TITLE: &str = "TD Emitter Locator";

pub const ROW_SPACING: f32 = 22.;
pub const COL_SPACING: f32 = 30.;

#[derive(Clone, Copy)]
enum Page {
    Sequence,
    Enzymes,
    Primters,
    Features,
}

pub fn draw(state: &mut State, ctx: &Context) {
    egui::CentralPanel::default().show(ctx, |ui| {
        // todo: Consider adding status like tshark/Wifi connection, monitor mode enabled, GPS coords etc.

        egui::containers::ScrollArea::vertical().show(ui, |ui| {
            ui.heading(RichText::new("System status").color(Color32::WHITE));
        });
    });
}
