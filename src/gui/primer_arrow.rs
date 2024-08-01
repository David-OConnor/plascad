//! This module contains code related to drawing primer arrows in the sequence view.

use std::ops::Range;

use eframe::{
    egui::{pos2, Align2, Color32, Direction, FontFamily, FontId, Pos2, Shape, Stroke, Ui},
    epaint::PathShape,
};

use crate::{
    gui::seq_view::{NT_WIDTH_PX, TEXT_X_START},
    primer::{
        PrimerData, PrimerDirection,
        PrimerDirection::{Forward, Reverse},
    },
};

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
    const STROKE_WIDTH: f32 = 2.;

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
        Stroke::new(STROKE_WIDTH, color_arrow),
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
            Stroke::new(STROKE_WIDTH, color_arrow),
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

/// Todo: Abstract this so it works for A: Arbitrary row count. B: Features.
// pub fn find_rows(primer_data: &PrimerData) {
pub fn find_rows(
    direction: PrimerDirection,
    seq_range: Range<usize>,
    ui: &mut Ui,
    nt_chars_per_row: usize,
    seq_len: usize,
    seq_i_to_px_rel: impl Fn(usize) -> Pos2,
) -> ((Pos2, Pos2), Option<(Pos2, Pos2)>) {
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
            }
        }
    };

    (bounds_row_0, bounds_row_1)
}

/// Add primer arrows to the display.
pub fn draw_primers(
    primer_data: &[PrimerData],
    ui: &mut Ui,
    nt_chars_per_row: usize,
    seq_len: usize,
    seq_i_to_px_rel: impl Fn(usize) -> Pos2,
) -> Vec<Shape> {
    let mut shapes = Vec::new();

    for prim_data in primer_data {
        let primer_matches = &prim_data.matches_seq;
        // todo: Sort out the direction. By matches, most likely.

        // todo: Do not run these calcs each time! Cache.
        for (direction, seq_range) in primer_matches {
            // todo: We need to match an arbitrary number of rows.
            let (bounds_row_0, bounds_row_1) = find_rows(
                *direction,
                seq_range.clone(),
                ui,
                nt_chars_per_row,
                seq_len,
                &seq_i_to_px_rel,
            );

            shapes.append(&mut primer_arrow(
                bounds_row_0,
                bounds_row_1,
                *direction,
                &prim_data.description,
                ui,
            ));
        }
    }
    shapes
}
