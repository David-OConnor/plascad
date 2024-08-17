//! This module contains GUI code related to the sequence view.

use std::process::id;

use eframe::egui::{
    text::CursorRange, Color32, Frame, RichText, ScrollArea, TextBuffer, TextEdit, Ui,
};

// todo: monospace font for all seqs.
use crate::sequence::{seq_from_str, seq_to_str, Feature, MIN_SEARCH_LEN};
use crate::{
    gui::{
        circle::feature_range_sliders,
        feature_table::feature_table,
        navigation::{page_seq_selector, page_seq_top_selector, PageSeq, PageSeqTop},
        primer_qc::primer_details,
        seq_view::sequence_vis,
        SPLIT_SCREEN_MAX_HEIGHT,
    },
    Selection,
};
// todo: monospace font for all seqs.
use crate::{
    gui::{COL_SPACING, ROW_SPACING},
    State,
};

fn seq_editor_raw(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        ui.heading("Sequence:");
        ui.label(&format!("len: {}", state.ui.seq_input.len()));
    });

    ScrollArea::vertical().id_source(200).show(ui, |ui| {
        let response = ui.add(TextEdit::multiline(&mut state.ui.seq_input).desired_width(800.));
        if response.changed() {
            state.generic.seq = seq_from_str(&state.ui.seq_input);
            state.ui.seq_input = seq_to_str(&state.generic.seq);
            state.sync_seq_related(None);
        }
    });
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
            ui.label(feature.location_descrip());
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
    ui.horizontal(|ui| {
        page_seq_top_selector(state, ui);

        // Sliders to edit the feature.
        feature_range_sliders(state, ui);
    });

    // Limit the top section height.
    let screen_height = ui.ctx().available_rect().height();
    let half_screen_height = screen_height / SPLIT_SCREEN_MAX_HEIGHT;

    match state.ui.page_seq_top {
        PageSeqTop::Primers => {
            Frame::none().show(ui, |ui| {
                ScrollArea::vertical()
                    .max_height(half_screen_height)
                    .id_source(69)
                    .show(ui, |ui| primer_details(state, ui));
            });
        }
        PageSeqTop::Features => {
            Frame::none().show(ui, |ui| {
                ScrollArea::vertical()
                    .max_height(half_screen_height)
                    .id_source(70)
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

        ui.label("🔍").on_hover_text(
            "Search the sequence and its complement for this term. (Ctrl + F to highlight)",
        );

        // This nonstandard way of adding the text input is required for the auto-highlight on ctrl+F behavior.
        let mut output = TextEdit::singleline(&mut state.ui.search_input)
            .desired_width(400.)
            .show(ui);
        let response = output.response;

        if state.ui.highlight_search_input {
            state.ui.highlight_search_input = false;
            response.request_focus();

            output.cursor_range = Some(CursorRange::select_all(&output.galley));
            // todo: Not working
        }

        if response.gained_focus() {
            state.ui.search_active = true; // Disable character entries in the sequence.
        }

        if response.changed {
            state.search_seq = seq_from_str(&state.ui.search_input);
            state.ui.search_input = seq_to_str(&state.search_seq); // Ensures only valid NTs are present.

            // todo: This still adds a single char, then blanks the cursor...
            state.ui.text_cursor_i = None; // Make sure we are not adding chars.
            state.sync_search();
        };

        if state.ui.search_input.len() >= MIN_SEARCH_LEN {
            let len = state.volatile.search_matches.len();
            let text = if len == 1 {
                "1 match".to_string()
            } else {
                format!("{} matches", len)
            };
            ui.label(text);
        }

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

    ScrollArea::vertical()
        .id_source(100)
        .show(ui, |ui| match state.ui.page_seq {
            PageSeq::EditRaw => {
                state.ui.text_cursor_i = None; // prevents double-edits
                seq_editor_raw(state, ui);
            }
            PageSeq::View => {
                sequence_vis(state, ui);
            }
        });
}
