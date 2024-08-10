//! This module contains GUI code related to the sequence view.

use eframe::egui::{Color32, RichText, TextEdit, Ui};

// todo: monospace font for all seqs.
use crate::sequence::{seq_from_str, seq_to_str, Feature};
use crate::{gui::{
    features::feature_table,
    navigation::{page_seq_selector, page_seq_top_selector, PageSeq, PageSeqTop},
    primer_qc::primer_details,
    seq_view::sequence_vis,
}, primer::make_cloning_primers, Selection};
// todo: monospace font for all seqs.
use crate::{
    gui::{COL_SPACING, ROW_SPACING},
    State,
};

fn seq_editor(state: &mut State, ui: &mut Ui) {
    // ui.heading("Amplification");

    // ui.add_space(ROW_SPACING);

    ui.horizontal(|ui| {
        ui.heading("Sequence:");
        ui.label(&format!("len: {}", state.ui.seq_input.len()));
    });

    let response = ui.add(TextEdit::multiline(&mut state.ui.seq_input).desired_width(800.));
    if response.changed() {
        state.generic.seq = seq_from_str(&state.ui.seq_input);
        state.ui.seq_input = seq_to_str(&state.generic.seq);
        state.sync_seq_related(None);
    }
}

fn seq_editor_slic(state: &mut State, ui: &mut Ui) {
    ui.heading("SLIC and FastCloning");

    ui.add_space(ROW_SPACING);

    ui.horizontal(|ui| {
        ui.label("Insert location: ");
        let mut entry = state.insert_loc.to_string();
        if ui
            .add(TextEdit::singleline(&mut entry).desired_width(40.))
            .changed()
        {
            state.insert_loc = entry.parse().unwrap_or(0);
            state.sync_cloning_product();
        }

        ui.add_space(COL_SPACING);

        if ui.button("➕ Make cloning primers").clicked() {
            make_cloning_primers(state);
        }

        if ui.button("Update seq with insert and vec").clicked() {
            state.sync_cloning_product();
            state.sync_seq_related(None);
        }
    });

    ui.horizontal(|ui| {
        ui.heading("Insert:");
        ui.label(&format!("len: {}", state.ui.seq_insert_input.len()));
    });

    let response = ui.add(
        TextEdit::multiline(&mut state.ui.seq_insert_input).desired_width(ui.available_width()),
    );
    if response.changed() {
        let seq = seq_from_str(&state.ui.seq_insert_input);
        state.ui.seq_insert_input = seq_to_str(&seq);
    }

    ui.add_space(ROW_SPACING);

    ui.horizontal(|ui| {
        ui.heading("Vector:");
        ui.label(&format!("len: {}", state.ui.seq_vector_input.len()));
    });

    let response = ui.add(
        TextEdit::multiline(&mut state.ui.seq_vector_input).desired_width(ui.available_width()),
    );
    if response.changed() {
        let seq = seq_from_str(&state.ui.seq_vector_input);
        state.ui.seq_vector_input = seq_to_str(&seq);
    }
}

/// Displays text of the feature under the cursor, or selected, as required.
fn feature_text(feature: &Option<usize>, features: &[Feature], ui: &mut Ui) {
    match feature {
        Some(i) => {
            if features.len() + 1 < *i {
                eprintln!("Invalid hover feature");
            }
            let feature = &features[*i];

            ui.label(&feature.label); // todo: IDeally heading here, but it is currnetly causing display jumping.
            ui.label(format!(
                "{}..{}",
                feature.index_range.0, feature.index_range.1
            ));
            let (r, g, b) = feature.color();
            ui.label(
                RichText::new(feature.feature_type.to_string()).color(Color32::from_rgb(r, g, b)),
            );

            // todo?
            for note in &feature.notes {
                // ui.label(&format!("{}: {}", note.0, note.1));
            }
        }
        None => (),
    }
}

/// Component for the sequence page.
pub fn seq_page(state: &mut State, ui: &mut Ui) {
    page_seq_top_selector(state, ui);

    match state.ui.page_seq_top {
        PageSeqTop::Primers => primer_details(state, ui),
        PageSeqTop::Features => feature_table(state, ui),
        PageSeqTop::None => (),
    }

    ui.add_space(ROW_SPACING);

    ui.horizontal(|ui| {
        page_seq_selector(state, ui);
        ui.add_space(COL_SPACING);

        let mut feature_to_disp = None;
        if let Selection::Feature(i) = state.ui.selected_item {
            feature_to_disp = Some(i);
        } else if state.ui.feature_hover.is_some() {
            feature_to_disp = Some(state.ui.feature_hover.unwrap());
        }

        feature_text(&feature_to_disp, &state.generic.features, ui);
    });

    ui.add_space(ROW_SPACING / 2.);

    match state.ui.page_seq {
        PageSeq::EditSeq => {
            seq_editor(state, ui);
        }
        PageSeq::EditSlic => {
            seq_editor_slic(state, ui);
        }
        PageSeq::View => {
            ui.horizontal(|ui| {
                // todo: DRY with above

                // todo: Impl insert loc changing A/R. Likely in a specialty cloning mode.
                // ui.label("Insert location: ");
                // let mut entry = state.insert_loc.to_string();
                // if ui.button("⏴").clicked() {
                //     if state.insert_loc > 0 {
                //         state.insert_loc -= 1;
                //     }
                //     state.sync_cloning_product();
                // };
                //
                // let response = ui.add(TextEdit::singleline(&mut entry).desired_width(40.));
                // if response.changed() {
                //     state.insert_loc = entry.parse().unwrap_or(0);
                //     state.sync_cloning_product();
                // }

                // if ui.button("⏵").clicked() {
                //     if state.insert_loc + 1 < state.ui.seq_vector_input.len() {
                //         state.insert_loc += 1;
                //     }
                //     state.sync_cloning_product();
                // };
            });

            sequence_vis(state, ui);
        }
    }
}
