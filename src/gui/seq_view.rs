//! This module contains GUI code related to the sequence visulization.

use std::ops::{Range, RangeInclusive};

use eframe::{
    egui::{
        pos2, vec2, Align2, Color32, FontFamily, FontId, Frame, Pos2, Rect, ScrollArea, Sense,
        Shape, Stroke, Ui,
    },
    emath::RectTransform,
    epaint::PathStroke,
};

use crate::{
    gui::{
        feature_from_index,
        feature_overlay::{draw_features, draw_selection},
        get_cursor_text,
        navigation::page_button,
        primer_arrow, COL_SPACING, ROW_SPACING,
    },
    sequence::ReadingFrame,
    util::{get_row_ranges, pixel_to_seq_i, seq_i_to_pixel},
    State, StateUi,
};

// Pub for use in `util` functions.
pub const FONT_SIZE_SEQ: f32 = 14.;
pub const COLOR_SEQ: Color32 = Color32::LIGHT_BLUE;
// (0xAD, 0xD8, 0xE6)
pub const COLOR_SEQ_DIMMED: Color32 = Color32::from_rgb(140, 160, 165); // Eg dim when there are search results
pub const COLOR_CODING_REGION: Color32 = Color32::from_rgb(255, 0, 170);
pub const COLOR_RE: Color32 = Color32::LIGHT_RED;
pub const COLOR_CURSOR: Color32 = Color32::from_rgb(255, 255, 0);
pub const COLOR_SEARCH_RESULTS: Color32 = Color32::from_rgb(255, 255, 130);
pub const COLOR_SELECTED_NTS: Color32 = Color32::from_rgb(255, 60, 255);

const BACKGROUND_COLOR: Color32 = Color32::from_rgb(10, 20, 10);

pub const NT_WIDTH_PX: f32 = 8.; // todo: Automatic way? This is valid for monospace font, size 14.
pub const VIEW_AREA_PAD_LEFT: f32 = 60.; // Bigger to accomodate the index display.
pub const VIEW_AREA_PAD_RIGHT: f32 = 20.;
// pub const SEQ_ROW_SPACING_PX: f32 = 34.;
pub const SEQ_ROW_SPACING_PX: f32 = 40.;

pub const TEXT_X_START: f32 = VIEW_AREA_PAD_LEFT;
pub const TEXT_Y_START: f32 = TEXT_X_START;
const MAX_SEQ_AREA_HEIGHT: u16 = 300;

/// These aguments define the circle, and are used in many places in this module.
pub struct SeqViewData {
    pub seq_len: usize,
    pub row_ranges: Vec<RangeInclusive<usize>>,
    pub to_screen: RectTransform,
    pub from_screen: RectTransform,
    // /// This is `from_screen * center`. We store it here to cache.
    // pub center_rel: Pos2,
}

impl SeqViewData {
    pub fn seq_i_to_px_rel(&self, i: usize) -> Pos2 {
        self.to_screen * seq_i_to_pixel(i, &self.row_ranges)
    }
}

fn draw_re_sites(state: &State, data: &SeqViewData, ui: &mut Ui) -> Vec<Shape> {
    let mut result = Vec::new();

    for (i_match, re_match) in state.volatile.restriction_enzyme_sites.iter().enumerate() {
        if re_match.lib_index + 1 > state.restriction_enzyme_lib.len() {
            continue;
        }
        let re = &state.restriction_enzyme_lib[re_match.lib_index];
        let cut_i = re_match.seq_index + 1; // to display in the right place.
        let cut_pos = data.seq_i_to_px_rel(cut_i + re.cut_after as usize);

        let bottom = pos2(cut_pos.x, cut_pos.y + 20.);

        result.push(Shape::LineSegment {
            points: [cut_pos, bottom],
            stroke: PathStroke::new(2., COLOR_RE),
        });

        let label_text = &re.name;
        let mut label_pos = pos2(cut_pos.x + 2., cut_pos.y - 4.);

        // Move the label position left if there is a nearby RE site on the right.
        // todo: Not appearing to be working well.
        let mut neighbor_on_right = false;
        for (i_other, re_match_other) in state.volatile.restriction_enzyme_sites.iter().enumerate()
        {
            if i_other != i_match
                && re_match_other.seq_index > re_match.seq_index
                && re_match_other.seq_index - re_match.seq_index < 10
            {
                neighbor_on_right = true;
                break;
            }
        }
        if neighbor_on_right {
            // label_pos.x -= 300.;
        }

        // Alternate above and below for legibility.
        // This requires the RE site list to be sorted by seq index, which it currently is.
        if i_match % 2 == 0 {
            label_pos.y += 30.;
        }

        // Add the label
        let label = ui.ctx().fonts(|fonts| {
            Shape::text(
                fonts,
                label_pos,
                Align2::LEFT_CENTER,
                label_text,
                FontId::new(16., FontFamily::Proportional),
                COLOR_RE,
            )
        });
        result.push(label)
    }
    result
}

/// Checkboxes to show or hide features.
pub fn display_filters(state_ui: &mut StateUi, ui: &mut Ui) {
    ui.horizontal(|ui| {
        ui.label("RE sites:");
        ui.checkbox(&mut state_ui.seq_visibility.show_res, "");
        ui.add_space(COL_SPACING / 2.);

        ui.label("Features:");
        ui.checkbox(&mut state_ui.seq_visibility.show_features, "");
        ui.add_space(COL_SPACING / 2.);

        ui.label("Primers:");
        ui.checkbox(&mut state_ui.seq_visibility.show_primers, "");
        ui.add_space(COL_SPACING / 2.);

        ui.label("Reading frame:");
        ui.checkbox(&mut state_ui.seq_visibility.show_reading_frame, "");
        ui.add_space(COL_SPACING / 2.);
    });
}

/// Draw each row's start sequence range to its left.
fn draw_seq_indexes(data: &SeqViewData, ui: &mut Ui) -> Vec<Shape> {
    let mut result = Vec::new();
    for range in &data.row_ranges {
        let mut pos = data.seq_i_to_px_rel(*range.start());
        pos.x -= VIEW_AREA_PAD_LEFT;

        // 1-based indexing.
        let text = range.start() + 1; // tood: 1k etc?

        result.push(ui.ctx().fonts(|fonts| {
            Shape::text(
                fonts,
                pos,
                Align2::LEFT_TOP,
                text,
                // Note: Monospace is important for sequences.
                FontId::new(FONT_SIZE_SEQ, FontFamily::Proportional),
                Color32::WHITE,
            )
        }));
    }

    result
}

fn orf_selector(state: &mut State, ui: &mut Ui) {
    ui.label("Reading frame:");

    let orf = &mut state.reading_frame;

    let orig = *orf;

    page_button(orf, ReadingFrame::Fwd0, ui, false);
    page_button(orf, ReadingFrame::Fwd1, ui, false);
    page_button(orf, ReadingFrame::Fwd2, ui, false);
    page_button(orf, ReadingFrame::Rev0, ui, false);
    page_button(orf, ReadingFrame::Rev1, ui, false);
    page_button(orf, ReadingFrame::Rev2, ui, false);

    if *orf != orig {
        state.sync_reading_frame()
    }
}

/// Find the sequence index under the cursor, if it is over the sequence.
/// 0-based indexing.
fn find_cursor_i(cursor_pos: Option<(f32, f32)>, data: &SeqViewData) -> Option<usize> {
    match cursor_pos {
        Some(p) => {
            // We've had issues where cursor above the seq would be treated as first row.
            let pos_relative = data.from_screen * pos2(p.0, p.1);

            if pos_relative.x > 0. && pos_relative.y > 0. {
                pixel_to_seq_i(pos_relative, &data.row_ranges)
            } else {
                None
            }
        }
        None => None,
    }
}

/// Draw the DNA sequence nucleotide by nucleotide. This allows fine control over color, and other things.
fn draw_nts(state: &State, data: &SeqViewData, ui: &mut Ui) -> Vec<Shape> {
    let mut result = Vec::new();

    for (i, nt) in state.generic.seq.iter().enumerate() {
        let pos = data.seq_i_to_px_rel(i);

        let letter_color = {
            let mut r = COLOR_SEQ;

            let mut highlighted = false;
            if state.ui.seq_visibility.show_reading_frame {
                for rf in &state.volatile.reading_frame_matches {
                    let (start, end) = rf.range;
                    let i_ = i + 1; // 1-based indexing
                    if start <= i_ && i_ <= end {
                        r = COLOR_CODING_REGION;
                        highlighted = true;
                    }
                }
            }

            // todo: We have inconsistencies in how we index across the board.
            // todo: Try resolving using the inclusive range type, and standardizing to 1-based.
            // todo: Resolve this one field at a time, from a working state.
            if let Some((start, end)) = state.ui.text_selection {
                let i_ = i + 1;
                if start <= i_ && i_ < end {
                    r = COLOR_SELECTED_NTS;
                }
            }

            // This overrides reading frame matches and text selection.
            for search_result in &state.volatile.search_matches {
                if search_result.range.contains(&i) {
                    r = COLOR_SEARCH_RESULTS;
                    highlighted = true;
                }
            }

            // Dim normal text if there are search results.
            if state.volatile.search_matches.len() > 0 && !highlighted {
                r = COLOR_SEQ_DIMMED;
            }

            r
        };

        result.push(ui.ctx().fonts(|fonts| {
            Shape::text(
                fonts,
                pos,
                Align2::LEFT_TOP,
                nt.as_str(),
                // Note: Monospace is important for sequences.
                FontId::new(FONT_SIZE_SEQ, FontFamily::Monospace),
                letter_color,
            )
        }));
    }

    result
}

fn draw_text_cursor(cursor_i: Option<usize>, data: &SeqViewData) -> Vec<Shape> {
    let mut result = Vec::new();

    if let Some(i) = cursor_i {
        let mut top = data.seq_i_to_px_rel(i);
        top.y -= 3.;
        let bottom = pos2(top.x, top.y + 23.);

        result.push(Shape::line_segment(
            [top, bottom],
            Stroke::new(2., COLOR_CURSOR),
        ));
    }

    result
}

/// Draw the sequence with primers, insertion points, and other data visible, A/R
pub fn sequence_vis(state: &mut State, ui: &mut Ui) {
    let mut shapes = vec![];

    let seq_len = state.generic.seq.len();

    state.ui.nt_chars_per_row = ((ui.available_width()
        - (VIEW_AREA_PAD_LEFT + VIEW_AREA_PAD_RIGHT))
        / NT_WIDTH_PX) as usize;
    let row_ranges = get_row_ranges(seq_len, state.ui.nt_chars_per_row);

    let mouse_posit_lbl = get_cursor_text(state.ui.cursor_seq_i, seq_len);
    let text_posit_lbl = get_cursor_text(state.ui.text_cursor_i, seq_len);

    ui.horizontal(|ui| {
        orf_selector(state, ui);
        ui.add_space(COL_SPACING);

        display_filters(&mut state.ui, ui);
        ui.add_space(COL_SPACING);

        ui.label("Cursor:");
        ui.heading(text_posit_lbl);

        ui.label("Mouse:");
        ui.heading(mouse_posit_lbl);
    });

    ScrollArea::vertical().show(ui, |ui| {
        Frame::canvas(ui.style())
            .fill(BACKGROUND_COLOR)
            .show(ui, |ui| {
                let (response, _painter) = {
                    // Estimate required height, based on seq len.
                    let total_seq_height = row_ranges.len() as f32 * SEQ_ROW_SPACING_PX + 60.;

                    let height = total_seq_height;

                    let desired_size = vec2(ui.available_width(), height);
                    // ui.allocate_painter(desired_size, Sense::click())
                    ui.allocate_painter(desired_size, Sense::click_and_drag())
                };

                let to_screen = RectTransform::from_to(
                    Rect::from_min_size(Pos2::ZERO, response.rect.size()),
                    response.rect,
                );

                let from_screen = to_screen.inverse();

                let data = SeqViewData {
                    seq_len,
                    row_ranges,
                    to_screen,
                    from_screen,
                };

                let prev_cursor_i = state.ui.cursor_seq_i;
                state.ui.cursor_seq_i = find_cursor_i(state.ui.cursor_pos, &data);

                if prev_cursor_i != state.ui.cursor_seq_i {
                    state.ui.feature_hover =
                        feature_from_index(&state.ui.cursor_seq_i, &state.generic.features);
                }

                // Removed: We select cursor position instead now.
                // select_feature(state, &from_screen);

                // todo: Move this into a function A/R.
                if state.ui.click_pending_handle {
                    // This is set up so that a click outside the text area won't reset the cursor.
                    if state.ui.cursor_seq_i.is_some() {
                        state.ui.text_cursor_i = state.ui.cursor_seq_i;
                        state.ui.search_active = false;
                        state.ui.text_selection = None;
                    }
                    state.ui.click_pending_handle = false;
                }

                shapes.extend(draw_seq_indexes(&data, ui));

                if state.ui.seq_visibility.show_primers {
                    shapes.append(&mut primer_arrow::draw_primers(
                        &state.generic.primers,
                        &data,
                        ui,
                    ));
                }

                if state.ui.seq_visibility.show_features {
                    shapes.append(&mut draw_features(&state.generic.features, &data, ui));
                }

                if state.ui.seq_visibility.show_res {
                    shapes.append(&mut draw_re_sites(state, &data, ui));
                }

                if let Some(selection) = &state.ui.text_selection {
                    shapes.append(&mut draw_selection(selection, &data, ui));
                }

                // Draw nucleotides arfter the selection, so it shows through the fill.
                shapes.append(&mut draw_nts(state, &data, ui));

                shapes.append(&mut draw_text_cursor(state.ui.text_cursor_i, &data));

                ui.painter().extend(shapes);
            });
    });
}
