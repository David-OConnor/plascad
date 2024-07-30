//! This module contains GUI code related to the sequence visulization.

use eframe::{
    egui::{
        pos2, vec2, Align2, Color32, FontFamily, FontId, Frame, Pos2, Rect, ScrollArea, Sense,
        Shape, Stroke, TextEdit, Ui,
    },
    emath::RectTransform,
    epaint::{PathShape, PathStroke},
};

use crate::{
    gui::{features::feature_add_disp, int_field, COL_SPACING, ROW_SPACING},
    primer::{
        PrimerDirection,
        PrimerDirection::{Forward, Reverse},
    },
    sequence::Feature,
    util::{get_row_ranges, pixel_to_seq_i, seq_i_to_pixel},
    State, StateUi,
};

// Pub for use in `util` functions.
pub const FONT_SIZE_SEQ: f32 = 14.;
pub const COLOR_SEQ: Color32 = Color32::LIGHT_BLUE;
pub const COLOR_INSERT: Color32 = Color32::from_rgb(255, 127, 39);
pub const COLOR_RE: Color32 = Color32::LIGHT_RED;

pub const NT_WIDTH_PX: f32 = 8.; // todo: Automatic way? This is valid for monospace font, size 14.
pub const VIEW_AREA_PAD: f32 = 40.;
pub const SEQ_ROW_SPACING_PX: f32 = 34.;

pub const TEXT_X_START: f32 = VIEW_AREA_PAD / 2.;
pub const TEXT_Y_START: f32 = TEXT_X_START;
const MAX_SEQ_AREA_HEIGHT: u16 = 300;

/// Make a visual arrow for a primer. For  use inside a Frame::canvas.
/// Note: This function is fragile, and was constructed partly by trial and error.
fn primer_arrow(
    mut bounds_r0: (Pos2, Pos2),
    mut bounds_r1: Option<(Pos2, Pos2)>, // Assumes no more than two rows.
    direction: PrimerDirection,
    label: &str,
    ui: &mut Ui,
) -> Vec<Shape> {
    let color_arrow = match direction {
        Forward => Color32::from_rgb(255, 0, 255),
        Reverse => Color32::LIGHT_YELLOW,
    };

    let color_label = Color32::LIGHT_GREEN;
    let arrow_width = 2.;

    const VERTICAL_OFFSET: f32 = 14.; // Number of pixels above the sequence text.
    const LABEL_OFFSET: f32 = 7.;
    const HEIGHT: f32 = 16.;
    const SLANT: f32 = 20.; // slant different, in pixels, for the arrow.

    const V_OFFSET_REV: f32 = 2. * VERTICAL_OFFSET;

    bounds_r0.0.y -= VERTICAL_OFFSET;
    bounds_r0.1.y -= VERTICAL_OFFSET;
    if let Some(b) = bounds_r1.as_mut() {
        b.0.y -= VERTICAL_OFFSET;
        b.1.y -= VERTICAL_OFFSET;
    }

    let mut result = Vec::new();

    let ctx = ui.ctx();

    let mut top_left = bounds_r0.0;

    // Slant only if single-line.
    let mut top_right = if bounds_r1.is_none() {
        pos2(bounds_r0.1.x - SLANT, bounds_r0.1.y)
    } else {
        pos2(bounds_r0.1.x, bounds_r0.1.y)
    };
    let mut bottom_left = pos2(bounds_r0.0.x, bounds_r0.0.y + HEIGHT);
    let mut bottom_right = pos2(bounds_r0.1.x, bounds_r0.1.y + HEIGHT);

    if direction == Reverse {
        std::mem::swap(&mut top_left, &mut top_right);
        std::mem::swap(&mut bottom_left, &mut bottom_right);

        top_left.y += V_OFFSET_REV;
        top_right.y += V_OFFSET_REV;
        bottom_left.y += V_OFFSET_REV;
        bottom_right.y += V_OFFSET_REV;
    }

    let points = vec![top_left, bottom_left, bottom_right, top_right];

    result.push(Shape::Path(PathShape::closed_line(
        points,
        Stroke::new(arrow_width, color_arrow),
    )));

    if let Some(b) = bounds_r1 {
        let mut top_left = b.0;
        let mut bottom_left = pos2(b.0.x, b.0.y + HEIGHT);
        let mut bottom_right = pos2(b.1.x, b.1.y + HEIGHT);
        let mut top_right = pos2(b.1.x - SLANT, b.1.y);

        // todo: DRY.
        if direction == Reverse {
            top_right.x += SLANT;
            bottom_left.x += SLANT;

            std::mem::swap(&mut top_left, &mut top_right);
            std::mem::swap(&mut bottom_left, &mut bottom_right);

            top_left.y += V_OFFSET_REV;
            top_right.y += V_OFFSET_REV;
            bottom_left.y += V_OFFSET_REV;
            bottom_right.y += V_OFFSET_REV;
        }

        let points = vec![top_left, bottom_left, bottom_right, top_right];

        result.push(Shape::Path(PathShape::closed_line(
            points,
            Stroke::new(arrow_width, color_arrow),
        )));
    }

    let label_start_x = match direction {
        Forward => bounds_r0.0.x,
        Reverse => bounds_r0.1.x,
    } + LABEL_OFFSET;

    let label_pos = match direction {
        Forward => pos2(label_start_x, bounds_r0.0.y + LABEL_OFFSET),
        Reverse => pos2(label_start_x, bounds_r0.0.y + LABEL_OFFSET + V_OFFSET_REV),
    };

    let label = ctx.fonts(|fonts| {
        Shape::text(
            fonts,
            label_pos,
            Align2::LEFT_CENTER,
            label,
            FontId::new(16., FontFamily::Proportional),
            color_label,
        )
    });

    result.push(label);
    result
}

fn features(state: &State, ui: &mut Ui, seq_i_to_px_rel: impl Fn(usize) -> Pos2) -> Vec<Shape> {
    let mut result = Vec::new();

    for feature in &state.features {
        // name, color, index_range
    }
    result
}

fn re_sites(state: &State, ui: &mut Ui, seq_i_to_px_rel: impl Fn(usize) -> Pos2) -> Vec<Shape> {
    let mut result = Vec::new();

    for (i_match, re_match) in state.restriction_enzyme_sites.iter().enumerate() {
        if re_match.lib_index + 1 > state.restriction_enzyme_lib.len() {
            continue;
        }
        let re = &state.restriction_enzyme_lib[re_match.lib_index];

        let cut_i = re_match.seq_index + 1; // to display in the right place.

        // let cut_i = match re_match.direction {
        //     Forward => re_match.seq_index + 1,
        //     // todo: Rev may be deprecated.
        //     Reverse => seq_len - re_match.seq_index,
        // };
        let cut_pos = seq_i_to_px_rel(cut_i + re.cut_after as usize);

        let bottom = pos2(cut_pos.x, cut_pos.y + 20.);

        result.push(Shape::LineSegment {
            points: [cut_pos, bottom],
            stroke: PathStroke::new(2., COLOR_RE),
        });

        // let label_text = format!("{} - {}", re.name, re_match.seq_index);
        let label_text = format!("{}", re.name);
        let mut label_pos = pos2(cut_pos.x + 2., cut_pos.y - 4.);

        // Move the label position left if there is a nearby RE site on the right.
        // todo: Not appearing to be working well.
        let mut neighbor_on_right = false;
        for (i_other, re_match_other) in state.restriction_enzyme_sites.iter().enumerate() {
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
fn display_filters(state_ui: &mut StateUi, ui: &mut Ui) {
    ui.horizontal(|ui| {
        ui.label("Show: ");
        ui.add_space(COL_SPACING);

        ui.label("RE sites:");
        ui.checkbox(&mut state_ui.show_res, "");
        ui.add_space(COL_SPACING / 2.);

        ui.label("Primers:");
        ui.checkbox(&mut state_ui.show_primers, "");
        ui.add_space(COL_SPACING / 2.);
    });
}

/// Add primer arrows to the display.
fn draw_primers(
    state: &State,
    ui: &mut Ui,
    nt_chars_per_row: usize,
    seq_len: usize,
    seq_i_to_px_rel: impl Fn(usize) -> Pos2,
) -> Vec<Shape> {
    let mut shapes = Vec::new();

    for prim_data in &state.primer_data {
        // let primer_matches = match state.ui.page_primer {
        //     PagePrimer::Amplification => &prim_data.matches_seq,
        //     PagePrimer::SlicFc => &prim_data.matches_vector_with_insert,
        //     PagePrimer::SlicFc => &prim_data.matches_vector_with_insert,
        // };
        let primer_matches = &prim_data.matches_seq;
        // todo: Sort out the direction. By matches, most likely.

        // todo: Do not run these calcs each time! Cache.
        for (direction, seq_range) in primer_matches {
            let (start, end) = match direction {
                Forward => (seq_range.start, seq_range.end),
                Reverse => (seq_len - seq_range.start, seq_len - seq_range.end),
            };

            let start_pos = seq_i_to_px_rel(start);
            let end_pos = seq_i_to_px_rel(end);

            // Check if we split across rows.
            let (bounds_row_0, bounds_row_1) = if start_pos.y == end_pos.y {
                ((start_pos, end_pos), None)
            } else {
                // let (col, row) = seq_i_to_col_row(seq_range.start);

                // let row_0_end = seq_i_to_pixel_rel(seq_range.start);

                match direction {
                    Forward => {
                        let row_0_end = pos2(
                            TEXT_X_START + NT_WIDTH_PX * (1. + nt_chars_per_row as f32),
                            start_pos.y,
                        );
                        let row_1_start = pos2(TEXT_X_START, end_pos.y);

                        ((start_pos, row_0_end), Some((row_1_start, end_pos)))
                    }
                    Reverse => {
                        // todo: DRY
                        // let row_0_end = pos2(
                        //     TEXT_X_START + NT_WIDTH_PX * (1. + nt_chars_per_row as f32),
                        //     start_pos.y,
                        // );

                        let row_0_end = pos2(
                            ui.available_width()
                                - (TEXT_X_START + NT_WIDTH_PX * (1. + nt_chars_per_row as f32)),
                            start_pos.y,
                        );

                        let row_1_start = pos2(ui.available_width(), end_pos.y);

                        ((row_0_end, start_pos), Some((end_pos, row_1_start)))
                        // ((start_pos, row_0_end), Some((row_1_start, end_pos)))
                    }
                }
            };

            let mut arrow_fwd = primer_arrow(
                bounds_row_0,
                bounds_row_1,
                *direction,
                &prim_data.description,
                ui,
            );
            shapes.append(&mut arrow_fwd);
        }
    }
    shapes
}

/// Draw the sequence with primers, insertion points, and other data visible, A/R
pub fn sequence_vis(state: &mut State, ui: &mut Ui) {
    let mut shapes = vec![];

    let seq_len = state.seq.len();

    let nt_chars_per_row = ((ui.available_width() - VIEW_AREA_PAD) / NT_WIDTH_PX) as usize; // todo: +1 etc?
    let row_ranges = get_row_ranges(seq_len, nt_chars_per_row);

    let cursor_posit_text = match state.ui.cursor_seq_i {
        Some(p) => {
            if p + 1 <= state.seq.len() {
                // + 1, as the convention is to use 1-based indexing vice 0.
                &(p + 1).to_string()
                // This occurs if the cursor is on the last row, right of the last NT.
            } else {
                ""
            }
        }
        None => "",
    };

    ui.horizontal(|ui| {
        feature_add_disp(state, ui);

        ui.add_space(COL_SPACING * 3.);

        ui.label(format!("Cursor: {}", cursor_posit_text));
    });

    display_filters(&mut state.ui, ui);

    ScrollArea::vertical().id_source(0).show(ui, |ui| {
        Frame::canvas(ui.style())
            .fill(Color32::from_rgb(10, 20, 10))
            .show(ui, |ui| {
                let (mut response, painter) = {
                    // Estimate required height, based on seq len.
                    let total_seq_height = row_ranges.len() as f32 * SEQ_ROW_SPACING_PX + 60.;
                    // leto height = min(total_seq_height as u16, MAX_SEQ_AREA_HEIGHT);

                    let height = total_seq_height;

                    let desired_size = vec2(ui.available_width(), height);
                    ui.allocate_painter(desired_size, Sense::click())
                };

                // let to_screen =
                //     RectTransform::from_to(Rect::from_x_y_ranges(0.0..=1.0, -1.0..=1.0), rect);

                let to_screen = RectTransform::from_to(
                    Rect::from_min_size(Pos2::ZERO, response.rect.size()),
                    response.rect,
                );

                let from_screen = to_screen.inverse();

                let seq_i_to_px_rel = |i| to_screen * seq_i_to_pixel(i, &row_ranges);

                // todo: Is this an apt place to handle this?
                {
                    state.ui.cursor_seq_i = match state.ui.cursor_pos {
                        Some(p) => {
                            // We've had issues where cursor above the seq would be treated as first row.
                            let pos_relative = from_screen * pos2(p.0, p.1);

                            if pos_relative.x > 0. && pos_relative.y > 0. {
                                pixel_to_seq_i(pos_relative, &row_ranges).map(|v| v)
                            } else {
                                None
                            }
                        }
                        None => None,
                    };
                }

                let ctx = ui.ctx();

                // Draw the sequence NT by NT. This allows fine control over color, and other things.
                for (i, nt) in state.seq.iter().enumerate() {
                    let pos = seq_i_to_px_rel(i);

                    let mut color = COLOR_SEQ;
                    // match state.ui.page_primer {
                    //     PagePrimer::Amplification => {}
                    //     PagePrimer::SlicFc => {
                    //         if i < state.insert_loc {
                    //         } else if i < state.insert_loc + state.ui.seq_insert_input.len() {
                    //             color = COLOR_INSERT;
                    //         } else {
                    //         }
                    //     }
                    // }

                    // todo: Needs an update
                    // if i < state.insert_loc + state.ui.seq_insert_input.len() {
                    //     color = COLOR_INSERT;
                    // }

                    shapes.push(ctx.fonts(|fonts| {
                        Shape::text(
                            fonts,
                            pos,
                            Align2::LEFT_TOP,
                            nt.as_str(),
                            // Note: Monospace is important for sequences.
                            FontId::new(FONT_SIZE_SEQ, FontFamily::Monospace),
                            color,
                        )
                    }));
                }

                if state.ui.show_primers {
                    shapes.append(&mut draw_primers(
                        state,
                        ui,
                        nt_chars_per_row,
                        seq_len,
                        seq_i_to_px_rel,
                    ));
                }

                if state.ui.show_res {
                    shapes.append(&mut re_sites(state, ui, seq_i_to_px_rel));
                }
                shapes.append(&mut features(state, ui, seq_i_to_px_rel));
                ui.painter().extend(shapes);
            });
    });

    ui.add_space(ROW_SPACING);
}
