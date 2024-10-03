use eframe::egui::{Color32, RichText, TextEdit, Ui};

use crate::{
    alignment::{align_pairwise, distance, DistanceType},
    gui::{theme::COLOR_ACTION, ROW_SPACING},
    sequence::{seq_from_str, seq_to_str, Nucleotide},
    State,
};

pub fn alignment_page(state: &mut State, ui: &mut Ui) {
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
        state.alignment.alignment_result = Some(align_pairwise(
            &state.alignment.seq_a,
            &state.alignment.seq_b,
        ));
        state.alignment.dist_result = Some(distance(
            &state.alignment.seq_a,
            &state.alignment.seq_b,
            DistanceType::Hamming,
        ));
    }

    ui.add_space(ROW_SPACING);
    if let Some(alignment) = &state.alignment.alignment_result {
        ui.heading(&format!("Alignment score: {:?}", alignment.score));
    }
    if let Some(dist) = &state.alignment.dist_result {
        ui.heading(&format!("Distance. Score: {:?}", dist));
    }
}
