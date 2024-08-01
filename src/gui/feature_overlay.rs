//! This module is related to drawing features on the sequence view. It is similar to `primer_arrow.rs`.

// todo: Abstract out diffs between this and the primer arrow; avoid repeated code.

use eframe::{
    egui::{pos2, Align2, Color32, FontFamily, FontId, Pos2, Shape, Stroke, Ui},
    epaint::PathShape,
};

use crate::sequence::{
    Feature,
    FeatureDirection::{Forward, Reverse},
    FeatureType,
};

/// Draw an indicator near a feature's sequence highlighting the feature visually, showing its label etc.
pub fn feature_overlay(
    feature: &Feature,
    mut bounds_r0: (Pos2, Pos2),
    mut bounds_r1: Option<(Pos2, Pos2)>, // Assumes no more than two rows.
    ui: &mut Ui,
) -> Vec<Shape> {
    // todo: Currently a C+P from primer arrow.

    let color_arrow = match feature.feature_type {
        FeatureType::Generic => Color32::from_rgb(255, 0, 255),
        FeatureType::Gene => Color32::from_rgb(255, 128, 128),
        FeatureType::Ori => Color32::from_rgb(255, 0, 255),
        FeatureType::RnaPolyBindSite => Color32::from_rgb(255, 0, 255),
    };
    // todo: Consider how you want to handle colors. Ie, perhaps defauilt to type colors, but allow
    // todo overriding.

    let color_label = Color32::LIGHT_GREEN;
    let arrow_width = 2.;

    const VERTICAL_OFFSET: f32 = 14.; // Number of pixels above the sequence text.
    const LABEL_OFFSET: f32 = 7.;
    const HEIGHT: f32 = 16.;

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
    let mut top_right = pos2(bounds_r0.1.x, bounds_r0.1.y);
    let mut bottom_left = pos2(bounds_r0.0.x, bounds_r0.0.y + HEIGHT);
    let mut bottom_right = pos2(bounds_r0.1.x, bounds_r0.1.y + HEIGHT);

    if feature.direction == Reverse {
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
        // let mut top_right = pos2(b.1.x - crate::gui::seq_view::SLANT, b.1.y);
        let mut top_right = pos2(b.1.x, b.1.y);

        // todo: DRY.
        // if direction == Reverse {
        //     top_right.x += crate::gui::seq_view::SLANT;
        //     bottom_left.x += crate::gui::seq_view::SLANT;
        //
        //     std::mem::swap(&mut top_left, &mut top_right);
        //     std::mem::swap(&mut bottom_left, &mut bottom_right);
        //
        //     top_left.y += V_OFFSET_REV;
        //     top_right.y += V_OFFSET_REV;
        //     bottom_left.y += V_OFFSET_REV;
        //     bottom_right.y += V_OFFSET_REV;
        // }

        let points = vec![top_left, bottom_left, bottom_right, top_right];

        result.push(Shape::Path(PathShape::closed_line(
            points,
            Stroke::new(arrow_width, color_arrow),
        )));
    }

    // let label_start_x = match feature.direction {
    //     Forward => bounds_r0.0.x,
    //     Reverse => bounds_r0.1.x,
    //     None => bounds_r0.0.x, // todo ?
    // } + LABEL_OFFSET;
    let label_start_x = 0.; // todo temp

    // let label_pos = match feature.direction {
    //     Forward => pos2(label_start_x, bounds_r0.0.y + LABEL_OFFSET),
    //     Reverse => pos2(label_start_x, bounds_r0.0.y + LABEL_OFFSET + V_OFFSET_REV),
    //     None => pos2(label_start_x, bounds_r0.0.y + LABEL_OFFSET), // todo?
    // };
    let label_pos = pos2(0., 0.);

    let label = ctx.fonts(|fonts| {
        Shape::text(
            fonts,
            label_pos,
            Align2::LEFT_CENTER,
            &feature.label,
            FontId::new(16., FontFamily::Proportional),
            color_label,
        )
    });

    result.push(label);
    result
}

pub fn draw_features(
    features: &[Feature],
    ui: &mut Ui,
    nt_chars_per_row: usize,
    seq_len: usize,
    seq_i_to_px_rel: impl Fn(usize) -> Pos2,
) -> Vec<Shape> {
    let mut shapes = Vec::new();

    // for feature in features {
    //     // let primer_matches = match state.ui.page_primer {
    //     //     PagePrimer::Amplification => &prim_data.matches_seq,
    //     //     PagePrimer::SlicFc => &prim_data.matches_vector_with_insert,
    //     //     PagePrimer::SlicFc => &prim_data.matches_vector_with_insert,
    //     // };
    //     // todo: Sort out the direction. By matches, most likely.
    //
    //     // todo: Do not run these calcs each time! Cache.
    //     let (start, end) = match feature.direction {
    //         Forward => (seq_range.start, seq_range.end),
    //         Reverse => (seq_len - seq_range.start, seq_len - seq_range.end),
    //     };
    //
    //     let start_pos = seq_i_to_px_rel(start);
    //     let end_pos = seq_i_to_px_rel(end);
    //
    //     // Check if we split across rows.
    //     let (bounds_row_0, bounds_row_1) = if start_pos.y == end_pos.y {
    //         ((start_pos, end_pos), None)
    //     } else {
    //         // let (col, row) = seq_i_to_col_row(seq_range.start);
    //
    //         // let row_0_end = seq_i_to_pixel_rel(seq_range.start);
    //
    //         match direction {
    //             Forward => {
    //                 let row_0_end = pos2(
    //                     TEXT_X_START + NT_WIDTH_PX * (1. + nt_chars_per_row as f32),
    //                     start_pos.y,
    //                 );
    //                 let row_1_start = pos2(TEXT_X_START, end_pos.y);
    //
    //                 ((start_pos, row_0_end), Some((row_1_start, end_pos)))
    //             }
    //             Reverse => {
    //                 // todo: DRY
    //                 // let row_0_end = pos2(
    //                 //     TEXT_X_START + NT_WIDTH_PX * (1. + nt_chars_per_row as f32),
    //                 //     start_pos.y,
    //                 // );
    //
    //                 let row_0_end = pos2(
    //                     ui.available_width()
    //                         - (TEXT_X_START + NT_WIDTH_PX * (1. + nt_chars_per_row as f32)),
    //                     start_pos.y,
    //                 );
    //
    //                 let row_1_start = pos2(ui.available_width(), end_pos.y);
    //
    //                 ((row_0_end, start_pos), Some((end_pos, row_1_start)))
    //                 // ((start_pos, row_0_end), Some((row_1_start, end_pos)))
    //             }
    //         }
    //     };

    //     shapes.append(&mut feature_overlay(
    //         feature,
    //         bounds_row_0,
    //         bounds_row_1,
    //         *direction,
    //         ui,
    //     ));
    // }
    shapes
}
