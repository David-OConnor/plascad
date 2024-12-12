use eframe::egui::{Color32, FontFamily, FontId, RichText, ScrollArea, TextEdit, Ui};
use na_seq::{seq_aa_from_str, seq_aa_to_str, seq_from_str, seq_to_str_lower};

use crate::{
    alignment::{align_pairwise_aa, align_pairwise_nt, distance_aa, distance_nt, AlignmentMode},
    gui::{
        theme::{COLOR_ACTION, COLOR_INFO},
        COL_SPACING, ROW_SPACING,
    },
    state::State,
};

fn mode_btn(state: &mut State, mode: AlignmentMode, name: &str, ui: &mut Ui) {
    let color = if state.alignment.mode == mode {
        Color32::LIGHT_BLUE
    } else {
        Color32::WHITE
    };
    if ui.button(RichText::new(name).color(color)).clicked() {
        state.alignment.mode = mode;

        // Reset fields to prevent invalid data.
        state.alignment.seq_a = Vec::new();
        state.alignment.seq_aa_a = Vec::new();
        state.alignment.seq_a_input = String::new();
        state.alignment.seq_b = Vec::new();
        state.alignment.seq_aa_b = Vec::new();
        state.alignment.seq_b_input = String::new();
    }
}

fn input_area(state: &mut State, seq_b: bool, ui: &mut Ui) {
    let (seq, seq_aa, seq_input) = if seq_b {
        (
            &mut state.alignment.seq_a,
            &mut state.alignment.seq_aa_a,
            &mut state.alignment.seq_a_input,
        )
    } else {
        (
            &mut state.alignment.seq_b,
            &mut state.alignment.seq_aa_b,
            &mut state.alignment.seq_b_input,
        )
    };

    let active_seq = state.generic[state.active].seq.clone(); // Avoid borrowing state in the closure

    ui.horizontal(|ui| {
        ui.label(if seq_b { "Seq B" } else { "Seq A" });
        ui.add_space(COL_SPACING);

        if let AlignmentMode::Dna = state.alignment.mode {
            if ui
                .button(RichText::new("Load active tab's sequence").color(Color32::LIGHT_BLUE))
                .clicked()
            {
                *seq = active_seq.clone(); // Use the pre-fetched active sequence
                *seq_input = seq_to_str_lower(seq);
            }
        }
    });

    let response = ui.add(TextEdit::multiline(seq_input).desired_width(800.));
    if response.changed() {
        match state.alignment.mode {
            AlignmentMode::Dna => {
                *seq = seq_from_str(seq_input);
                *seq_input = seq_to_str_lower(seq);
                // Todo: why?
                // state.sync_seq_related(None);
            }
            AlignmentMode::AminoAcid => {
                *seq_aa = seq_aa_from_str(seq_input);
                *seq_input = seq_aa_to_str(seq_aa);
            }
        }
    }
}

pub fn alignment_page(state: &mut State, ui: &mut Ui) {
    ui.add_space(ROW_SPACING);

    ui.horizontal(|ui| {
        ui.heading("Alignment (Work in Progress)");

        ui.add_space(COL_SPACING * 2.);

        mode_btn(state, AlignmentMode::Dna, "DNA", ui);
        mode_btn(state, AlignmentMode::AminoAcid, "Amino acid", ui);
    });
    ui.add_space(ROW_SPACING);

    ScrollArea::vertical().id_salt(200).show(ui, |ui| {
        input_area(state, false, ui);
        ui.add_space(ROW_SPACING);
        input_area(state, true, ui);

        ui.add_space(ROW_SPACING);

        if ui
            .button(RichText::new("Align").color(COLOR_ACTION))
            .clicked()
        {
            let alignment = match state.alignment.mode {
                AlignmentMode::Dna => {
                    align_pairwise_nt(&state.alignment.seq_a, &state.alignment.seq_b)
                }
                AlignmentMode::AminoAcid => {
                    align_pairwise_aa(&state.alignment.seq_aa_a, &state.alignment.seq_aa_b)
                }
            };

            state.alignment.alignment_result = Some(alignment.0);
            state.alignment.text_display = alignment.1;

            state.alignment.dist_result = Some(match state.alignment.mode {
                AlignmentMode::Dna => distance_nt(&state.alignment.seq_a, &state.alignment.seq_b),
                AlignmentMode::AminoAcid => {
                    distance_aa(&state.alignment.seq_aa_a, &state.alignment.seq_aa_b)
                }
            });
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
    });
}
