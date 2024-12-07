use eframe::egui::{Color32, FontFamily, FontId, RichText, ScrollArea, TextEdit, Ui};
use na_seq::{seq_from_str, seq_to_str};

use crate::{alignment::{align_pairwise, distance}, AlignmentMode, gui::{
    theme::{COLOR_ACTION, COLOR_INFO},
    ROW_SPACING,
}, State};
use crate::gui::COL_SPACING;

pub fn alignment_page(state: &mut State, ui: &mut Ui) {

    ui.add_space(ROW_SPACING);

    ui.horizontal(|ui| {
        ui.heading("Alignment (Work in Progress)");

        ui.add_space(COL_SPACING * 2.);

        let color = if state.alignment.mode == AlignmentMode::Dna {
            Color32::LIGHT_BLUE
        } else {
            Color32::WHITE
        };
       if ui.button(RichText::new("DNA").color(color)).clicked() {
            state.alignment.mode = AlignmentMode::Dna;
       }

        let color = if state.alignment.mode == AlignmentMode::AminoAcid {
            Color32::LIGHT_BLUE
        } else {
            Color32::WHITE
        };
        if ui.button(RichText::new("Amino acid").color(color)).clicked() {
            state.alignment.mode = AlignmentMode::AminoAcid;
        }
    });
    ui.add_space(ROW_SPACING);

    ScrollArea::vertical().id_salt(200).show(ui, |ui| {
        ui.horizontal(|ui| {
            ui.label("Seq A:");
            ui.add_space(COL_SPACING);

            if ui.button(RichText::new("Load active tab's sequence").color(Color32::LIGHT_BLUE)).clicked() {
                state.alignment.seq_a = state.generic[state.active].seq.clone();
                state.alignment.seq_a_input = seq_to_str(&state.alignment.seq_a);
            }
        });

        let response =
            ui.add(TextEdit::multiline(&mut state.alignment.seq_a_input).desired_width(800.));
        if response.changed() {
            state.alignment.seq_a = seq_from_str(&state.alignment.seq_a_input);
            state.alignment.seq_a_input = seq_to_str(&state.alignment.seq_a);
            state.sync_seq_related(None);
        }
        ui.add_space(ROW_SPACING);

        // todo: DRY. Also with sequence.mod.
        ui.horizontal(|ui| {
            ui.label("Seq B:");
            ui.add_space(COL_SPACING);

            if ui.button(RichText::new("Load active tab's sequence").color(Color32::LIGHT_BLUE)).clicked() {
                state.alignment.seq_b = state.generic[state.active].seq.clone();
                state.alignment.seq_b_input = seq_to_str(&state.alignment.seq_b);
            }
        });

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
    });
}
