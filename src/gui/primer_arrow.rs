//! This module contains code related to drawing primer arrows in the sequence view.

use std::ops::Range;

use eframe::{
    egui::{pos2, Align2, Color32, Direction, FontFamily, FontId, Pos2, Shape, Stroke, Ui},
    epaint::PathShape,
};
use eframe::egui::accesskit::VerticalOffset;
use crate::{
    gui::seq_view::{NT_WIDTH_PX, TEXT_X_START},
    primer::{
        PrimerData, PrimerDirection,
        PrimerDirection::{Forward, Reverse},
    },
};

const STROKE_WIDTH: f32 = 2.;

const VERTICAL_OFFSET_PRIMER: f32 = 14.; // Number of pixels above the sequence text.
const LABEL_OFFSET: f32 = 7.;
const HEIGHT: f32 = 16.;
const SLANT: f32 = 20.; // slant different, in pixels, for the arrow.

/// Make a visual arrow for a primer. For  use inside a Frame::canvas.
/// Note: This function is fragile, and was constructed partly by trial and error.
pub fn primer_arrow(
    row_ranges_px: &[(Pos2, Pos2)],
    // mut bounds_r0: (Pos2, Pos2),
    // mut bounds_r1: Option<(Pos2, Pos2)>, // Assumes no more than two rows.
    vertical_offset: f32,
    direction: PrimerDirection,
    label: &str,
    ui: &mut Ui,
) -> Vec<Shape> {
    println!("RRP LEN: {}", row_ranges_px.len());

    if row_ranges_px.is_empty() {
        return Vec::new();
    }

    let color_arrow = match direction {
        Forward => Color32::from_rgb(255, 0, 255),
        Reverse => Color32::LIGHT_YELLOW,
    };

    let color_label = Color32::LIGHT_GREEN;

    let v_offset_rev = 2. * vertical_offset;

    // Apply a vertical offset from the sequence.
    // todo:
    let row_ranges_px: Vec<(Pos2, Pos2)> = row_ranges_px
        .iter()
        .map(|(start, end)| {
            (
                pos2(start.x, start.y - vertical_offset),
                pos2(end.x, end.y - vertical_offset),

            )
        })
        .collect();

    let mut result = Vec::new();

    println!("\n A");
    for (start, end) in &row_ranges_px {
        let mut top_left = *start;
        let mut top_right = *end;
        let bottom_left = pos2(start.x, start.y + HEIGHT);
        let bottom_right = pos2(end.x, end.y + HEIGHT);


        println!("Points: {:?}", [top_left, bottom_left, bottom_right, top_right]);
        // todo: Update with slant.

        result.push(Shape::Path(PathShape::closed_line(
            vec![top_left, bottom_left, bottom_right, top_right],
            Stroke::new(STROKE_WIDTH, color_arrow),
        )));
    }
    //
    //
    //
    // // Slant only if single-line.
    // let mut top_right = if bounds_r1.is_none() {
    //     pos2(bounds_r0.1.x - SLANT, bounds_r0.1.y)
    // } else {
    //     pos2(bounds_r0.1.x, bounds_r0.1.y)
    // };
    // let mut bottom_left = pos2(bounds_r0.0.x, bounds_r0.0.y + HEIGHT);
    // let mut bottom_right = pos2(bounds_r0.1.x, bounds_r0.1.y + HEIGHT);
    //
    // if direction == Reverse {
    //     std::mem::swap(&mut top_left, &mut top_right);
    //     std::mem::swap(&mut bottom_left, &mut bottom_right);
    //
    //     top_left.y += V_OFFSET_REV;
    //     top_right.y += V_OFFSET_REV;
    //     bottom_left.y += V_OFFSET_REV;
    //     bottom_right.y += V_OFFSET_REV;
    // }

    // if let Some(b) = bounds_r1 {
    //     let mut top_left = b.0;
    //     let mut bottom_left = pos2(b.0.x, b.0.y + HEIGHT);
    //     let mut bottom_right = pos2(b.1.x, b.1.y + HEIGHT);
    //     let mut top_right = pos2(b.1.x - SLANT, b.1.y);
    //
    //     // todo: DRY.
    //     if direction == Reverse {
    //         top_right.x += SLANT;
    //         bottom_left.x += SLANT;
    //
    //         std::mem::swap(&mut top_left, &mut top_right);
    //         std::mem::swap(&mut bottom_left, &mut bottom_right);
    //
    //         top_left.y += V_OFFSET_REV;
    //         top_right.y += V_OFFSET_REV;
    //         bottom_left.y += V_OFFSET_REV;
    //         bottom_right.y += V_OFFSET_REV;
    //     }
    //
    //     result.push(Shape::Path(PathShape::closed_line(
    //         vec![top_left, bottom_left, bottom_right, top_right],
    //         Stroke::new(STROKE_WIDTH, color_arrow),
    //     )));
    // }

    let label_start_x = match direction {
        Forward => row_ranges_px[0].0.x,
        Reverse => row_ranges_px[row_ranges_px.len() - 1].1.x,
    } + LABEL_OFFSET;

    let label_pos = match direction {
        Forward => pos2(label_start_x, row_ranges_px[0].0.y + LABEL_OFFSET),
        Reverse => pos2(
            label_start_x,
            row_ranges_px[0].0.y + LABEL_OFFSET + v_offset_rev,
        ),
    };

    let label = ui.ctx().fonts(|fonts| {
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


// todo: Abstract out/use in feature overlay
/// Given an index range of a feature, return sequence ranges for each row the feature occupies.
pub fn get_feature_ranges(range: &Range<usize>, row_ranges: &[Range<usize>]) -> Vec<Range<usize>> {

}

/// Add primer arrows to the display.
pub fn draw_primers(
    primer_data: &[PrimerData],
    row_ranges: &[Range<usize>],
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
            // let (bounds_row_0, bounds_row_1) = get_row_range_px(
            //     row_ranges,
            //     *direction,
            //     seq_range.clone(),
            //     ui,
            //     nt_chars_per_row,
            //     seq_len,
            //     &seq_i_to_px_rel,
            // );

            let feature_ranges = get_feature_ranges(seq_range, row_ranges);

            let row_ranges_px: Vec<(Pos2, Pos2)> = feature_ranges
                .iter()
                .map(|r| (seq_i_to_px_rel(r.start), seq_i_to_px_rel(r.end)))
                .collect();

            // todo: PUt back; temp check on compiling.
            shapes.append(&mut primer_arrow(
                &row_ranges_px,
                // bounds_row_0,
                // bounds_row_1,
                VERTICAL_OFFSET_PRIMER,
                *direction,
                &prim_data.description,
                ui,
            ));
        }
    }
    shapes
}
