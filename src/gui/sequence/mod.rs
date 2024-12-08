//! This module contains GUI code related to the sequence view.

use eframe::egui::{text::CursorRange, Color32, Frame, RichText, ScrollArea, TextEdit, Ui};
use na_seq::{seq_complement, seq_from_str, seq_to_str_lower};

// todo: monospace font for all seqs.
use crate::misc_types::{Feature, FeatureDirection, MIN_SEARCH_LEN};
use crate::{
    gui::{
        circle::feature_range_sliders,
        feature_table::{direction_picker, feature_table},
        navigation::{page_seq_selector, page_seq_top_selector, PageSeq, PageSeqTop},
        primer_table::primer_details,
        sequence::seq_view::sequence_vis,
        theme::COLOR_ACTION,
        PRIMER_FWD_COLOR, SPLIT_SCREEN_MAX_HEIGHT,
    },
    primer::{Primer, PrimerData},
    util::RangeIncl,
    Selection,
};
// todo: monospace font for all seqs.
use crate::{
    gui::{COL_SPACING, ROW_SPACING},
    State,
};

mod feature_overlay;
mod primer_overlay;
pub mod seq_view;

fn seq_editor_raw(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        ui.heading("Sequence:");
        ui.label(&format!("len: {}", state.ui.seq_input.len()));
    });

    ScrollArea::vertical().id_salt(200).show(ui, |ui| {
        let response = ui.add(TextEdit::multiline(&mut state.ui.seq_input).desired_width(800.));
        if response.changed() {
            state.generic[state.active].seq = seq_from_str(&state.ui.seq_input);
            state.ui.seq_input = seq_to_str_lower(state.get_seq());
            state.sync_seq_related(None);
        }
    });
}

/// Displays text of the feature under the cursor, or selected, as required.
fn feature_text(i: usize, features: &[Feature], seq_len: usize, ui: &mut Ui) {
    if i >= features.len() {
        eprintln!("Invalid selected feature");
        return; // todo: Ideally set the feature to none.
    }
    let feature = &features[i];

    ui.label(&feature.label);
    ui.label(feature.location_descrip(seq_len));
    let (r, g, b) = feature.color();
    ui.label(RichText::new(feature.feature_type.to_string()).color(Color32::from_rgb(r, g, b)));

    // todo?
    for note in &feature.notes {
        // ui.label(&format!("{}: {}", note.0, note.1));
    }
}

fn primer_text(i: usize, primers: &[Primer], seq_len: usize, ui: &mut Ui) {
    if i >= primers.len() {
        eprintln!("Invalid selected primer");
        return; // todo: Ideally set the feature to none.
    }
    let primer = &primers[i];

    ui.label(&primer.name);
    ui.label(&primer.location_descrip());
    // todo: Rev color A/R
    ui.label(RichText::new(seq_to_str_lower(&primer.sequence)).color(PRIMER_FWD_COLOR));

    ui.label(&primer.description.clone().unwrap_or_default());
}

/// Add a toolbar to create a feature from selection, if appropriate.
fn feature_from_sel(state: &mut State, ui: &mut Ui) {
    if let Some(text_sel) = state.ui.text_selection {
        if ui
            .button(RichText::new("➕ Add feature from sel").color(COLOR_ACTION))
            .clicked()
        {
            state.generic[state.active].features.push(Feature {
                range: text_sel,
                label: state.ui.quick_feature_add_name.clone(),
                direction: state.ui.quick_feature_add_dir,
                ..Default::default()
            });

            state.ui.text_selection = None;
            state.ui.quick_feature_add_name = String::new();
        }

        if ui
            .button(RichText::new("➕ Add primer from sel").color(COLOR_ACTION))
            .clicked()
        {
            // todo: DRY with genbank parsing; common fn A/R.
            let seq = state.get_seq();
            let compl = &seq_complement(seq);
            let seq_primer = match state.ui.quick_feature_add_dir {
                FeatureDirection::Reverse => {
                    let range = RangeIncl::new(
                        seq.len() - text_sel.end + 1,
                        seq.len() - text_sel.start + 1,
                    );

                    range.index_seq(compl).unwrap_or_default()
                }
                _ => text_sel.index_seq(seq).unwrap_or_default(),
            }
            .to_vec();

            let volatile = PrimerData::new(&seq_primer);

            state.generic[state.active].primers.push(Primer {
                sequence: seq_primer,
                name: state.ui.quick_feature_add_name.clone(),
                description: None,
                volatile,
            });

            state.ui.quick_feature_add_name = String::new();
            state.sync_primer_matches(None);
            state.sync_primer_metrics();
        }

        direction_picker(&mut state.ui.quick_feature_add_dir, 200, ui);

        ui.label("Name:");
        if ui
            .add(TextEdit::singleline(&mut state.ui.quick_feature_add_name).desired_width(80.))
            .gained_focus()
        {
            state.ui.text_edit_active = true; // Disable character entries in the sequence.
        }

        ui.add_space(COL_SPACING)
    }
}

/// Component for the sequence page.
pub fn seq_page(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        page_seq_top_selector(state, ui);

        ui.add_space(COL_SPACING);

        feature_from_sel(state, ui);

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
                    .id_salt(69)
                    .show(ui, |ui| primer_details(state, ui));
            });
        }
        PageSeqTop::Features => {
            Frame::none().show(ui, |ui| {
                ScrollArea::vertical()
                    .max_height(half_screen_height)
                    .id_salt(70)
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
            state.ui.text_edit_active = true; // Disable character entries in the sequence.
            response.request_focus();

            output.cursor_range = Some(CursorRange::select_all(&output.galley));
            // todo: Not working
        }

        if response.gained_focus() {
            state.ui.text_edit_active = true; // Disable character entries in the sequence.
            println!("GF");
        }

        if response.changed {
            state.ui.text_edit_active = true;
            state.search_seq = seq_from_str(&state.ui.search_input);
            state.ui.search_input = seq_to_str_lower(&state.search_seq); // Ensures only valid NTs are present.

            // todo: This still adds a single char, then blanks the cursor...
            state.ui.text_cursor_i = None; // Make sure we are not adding chars.
            state.sync_search();
        };

        if state.ui.search_input.len() >= MIN_SEARCH_LEN {
            let len = state.volatile[state.active].search_matches.len();
            let text = if len == 1 {
                "1 match".to_string()
            } else {
                format!("{} matches", len)
            };
            ui.label(text);
        }

        ui.add_space(COL_SPACING);

        let mut feature_to_disp = None;
        let mut primer_to_disp = None;

        match state.ui.selected_item {
            Selection::Feature(i) => feature_to_disp = Some(i),
            Selection::Primer(i) => primer_to_disp = Some(i),
            Selection::None => {
                if state.ui.feature_hover.is_some() {
                    feature_to_disp = Some(state.ui.feature_hover.unwrap());
                }
            }
        }

        if let Some(feature_i) = feature_to_disp {
            feature_text(
                feature_i,
                &state.generic[state.active].features,
                state.get_seq().len(),
                ui,
            );
        }

        if let Some(primer_i) = primer_to_disp {
            primer_text(
                primer_i,
                &state.generic[state.active].primers,
                state.get_seq().len(),
                ui,
            );
        }
    });

    ui.add_space(ROW_SPACING / 2.);

    ScrollArea::vertical()
        .id_salt(100)
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
