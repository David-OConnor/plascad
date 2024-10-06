use eframe::egui::{FontFamily, FontId, RichText, TextEdit, Ui};
use seq::{seq_from_str, seq_to_str};

use crate::{
    alignment::{align_pairwise, distance},
    gui::{
        theme::{COLOR_ACTION, COLOR_INFO},
        ROW_SPACING,
    },
    State,
};

pub fn alignment_page(state: &mut State, ui: &mut Ui) {
    ui.heading("Alignment (Work in Progress");
    ui.add_space(ROW_SPACING);

    ui.label("Seq A:");
    let response =
        ui.add(TextEdit::multiline(&mut state.alignment.seq_a_input).desired_width(800.));
    if response.changed() {
        state.alignment.seq_a = seq_from_str(&state.alignment.seq_a_input);
        state.alignment.seq_a_input = seq_to_str(&state.alignment.seq_a);
        state.sync_seq_related(None);
    }
    ui.add_space(ROW_SPACING);

    // todo: DRY. Also with sequence.mod.
    ui.label("Seq B:");
    let response =
        ui.add(TextEdit::multiline(&mut state.alignment.seq_b_input).desired_width(800.));
    if response.changed() {
        state.alignment.seq_b = seq_from_str(&state.alignment.seq_b_input);
        state.alignment.seq_b_input = seq_to_str(&state.alignment.seq_b);
        state.sync_seq_related(None);
    }
    ui.add_space(ROW_SPACING);

    if ui
        .button(RichText::new("Align").color(COLOR_ACTION))
        .clicked()
    {
        let alignment = align_pairwise(&state.alignment.seq_a, &state.alignment.seq_b);
        state.alignment.alignment_result = Some(alignment.0);
        state.alignment.text_display = alignment.1;

        state.alignment.dist_result =
            Some(distance(&state.alignment.seq_a, &state.alignment.seq_b));
    }

    ui.add_space(ROW_SPACING);
    if let Some(alignment) = &state.alignment.alignment_result {
        ui.heading(&format!("Alignment score: {:?}", alignment.score));

        ui.add_space(ROW_SPACING);

        ui.label(
            RichText::new(&state.alignment.text_display)
                .color(COLOR_INFO)
                .font(FontId::new(16., FontFamily::Monospace)),
        );
    }
    // ui.add_space(ROW_SPACING);

    if let Some(dist) = &state.alignment.dist_result {
        ui.horizontal(|ui| {
            ui.heading(&format!("Distance score: {:?}", dist));
            let dist_type_text = if state.alignment.seq_a.len() == state.alignment.seq_b.len() {
                "(Hamming)"
            } else {
                "(Levenshtein)"
            };
            ui.label(dist_type_text);
        });
    }
}
