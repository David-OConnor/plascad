//! This module contains code related to drawing primer arrows in the sequence view.

use std::ops::Range;

use eframe::{
    egui::{pos2, Align2, Color32, FontFamily, FontId, Pos2, Shape, Stroke, Ui},
    epaint::PathShape,
};

use crate::{
    gui::seq_view::NT_WIDTH_PX,
    primer::{PrimerData, PrimerDirection},
    sequence::FeatureDirection::{self, Forward, Reverse},
    util,
};

const STROKE_WIDTH: f32 = 2.;

const VERTICAL_OFFSET_PRIMER: f32 = 14.; // Number of pixels above the sequence text.
const LABEL_OFFSET: f32 = 7.;
const HEIGHT: f32 = 16.;
const SLANT: f32 = 20.; // slant different, in pixels, for the arrow.

/// Make a visual arrow for a primer. For  use inside a Frame::canvas.
/// Note: This function is fragile, and was constructed partly by trial and error.
/// todo: Reanme draw_feature_overlay etc.
pub fn primer_arrow(
    feature_ranges_px: &[(Pos2, Pos2)],
    vertical_offset: f32,
    direction: FeatureDirection,
    label: &str,
    ui: &mut Ui,
) -> Vec<Shape> {
    if feature_ranges_px.is_empty() {
        return Vec::new();
    }

    let outline_color = match direction {
        Forward => Color32::from_rgb(255, 0, 255),
        Reverse => Color32::LIGHT_YELLOW,
        FeatureDirection::None => Color32::GOLD,
    };

    let color_label = Color32::LIGHT_GREEN;

    let v_offset_rev = 2. * vertical_offset;

    // Apply a vertical offset from the sequence.
    let feature_ranges_px: Vec<(Pos2, Pos2)> = feature_ranges_px
        .iter()
        .map(|(start, end)| {
            (
                pos2(start.x, start.y - vertical_offset),
                pos2(end.x, end.y - vertical_offset),
            )
        })
        .collect();

    let mut result = Vec::new();

    for (start, end) in &feature_ranges_px {
        let mut top_left = *start;
        let mut top_right = pos2(end.x + NT_WIDTH_PX, end.y);
        let bottom_left = pos2(start.x, start.y + HEIGHT);
        let bottom_right = pos2(end.x + NT_WIDTH_PX, end.y + HEIGHT);

        // todo: Update with slant.

        result.push(Shape::Path(PathShape::closed_line(
            vec![top_left, bottom_left, bottom_right, top_right],
            Stroke::new(STROKE_WIDTH, outline_color),
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

    // todo: Examine.
    let label_start_x = match direction {
        Forward => feature_ranges_px[0].0.x,
        Reverse => feature_ranges_px[feature_ranges_px.len() - 1].1.x,
        FeatureDirection::None => feature_ranges_px[0].0.x,
    } + LABEL_OFFSET;

    let label_pos = match direction {
        Forward => pos2(label_start_x, feature_ranges_px[0].0.y + LABEL_OFFSET),
        Reverse => pos2(
            label_start_x,
            feature_ranges_px[0].0.y + LABEL_OFFSET + v_offset_rev,
        ),
        FeatureDirection::None => pos2(label_start_x, feature_ranges_px[0].0.y + LABEL_OFFSET), // todo: Examine
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

/// Add primer arrows to the display.
pub fn draw_primers(
    primer_data: &[PrimerData],
    row_ranges: &[Range<usize>],
    ui: &mut Ui,
    seq_i_to_px_rel: impl Fn(usize) -> Pos2,
) -> Vec<Shape> {
    let mut shapes = Vec::new();

    for prim_data in primer_data {
        let primer_matches = &prim_data.matches_seq;

        // todo: Do not run these calcs each time. Cache.
        for (direction, seq_range) in primer_matches {
            let feature_ranges = util::get_feature_ranges(seq_range, row_ranges);

            let feature_ranges_px: Vec<(Pos2, Pos2)> = feature_ranges
                .iter()
                .map(|r| (seq_i_to_px_rel(r.start), seq_i_to_px_rel(r.end)))
                .collect();

            // todo: PUt back; temp check on compiling.
            shapes.append(&mut primer_arrow(
                &feature_ranges_px,
                VERTICAL_OFFSET_PRIMER,
                (*direction).into(),
                &prim_data.description,
                ui,
            ));
        }
    }
    shapes
}
