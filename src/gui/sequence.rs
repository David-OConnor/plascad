//! This module contains GUI code related to the sequence view.

use eframe::egui::{Color32, Frame, RichText, ScrollArea, TextEdit, Ui};

// todo: monospace font for all seqs.
use crate::sequence::{seq_from_str, seq_to_str, Feature};
use crate::{
    gui::{
        SPLIT_SCREEN_MAX_HEIGHT,
        features::feature_table,
        navigation::{page_seq_selector, page_seq_top_selector, PageSeq, PageSeqTop},
        primer_qc::primer_details,
        seq_view::sequence_vis,
    },
    Selection,
};
// todo: monospace font for all seqs.
use crate::{
    gui::{COL_SPACING, ROW_SPACING},
    State,
};
fn seq_editor(state: &mut State, ui: &mut Ui) {
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

/// Displays text of the feature under the cursor, or selected, as required.
fn feature_text(feature: &Option<usize>, features: &[Feature], ui: &mut Ui) {
    match feature {
        Some(i) => {
            if features.len() < *i + 1 {
                eprintln!("Invalid hover feature");
                return; // todo: Ideally set the feature to none.
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

    // Limit the top section height.
    let screen_height = ui.ctx().available_rect().height();
    let half_screen_height = screen_height / SPLIT_SCREEN_MAX_HEIGHT;

    match state.ui.page_seq_top {
        PageSeqTop::Primers => {
            Frame::none().show(ui, |ui| {
                ScrollArea::vertical()
                    .max_height(half_screen_height)
                    .show(ui, |ui| primer_details(state, ui));
            });
        }
        PageSeqTop::Features => {
            Frame::none().show(ui, |ui| {
                ScrollArea::vertical()
                    .max_height(half_screen_height)
                    .show(ui, |ui| {
                        feature_table(state, ui);
                    });
            });
        }
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
        // PageSeq::EditSlic => {
        //     cloning::seq_editor_slic(state, ui);
        // }
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
