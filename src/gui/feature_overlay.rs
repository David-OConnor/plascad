//! This module is related to drawing features on the sequence view. It is similar to `primer_arrow.rs`.

// todo: Abstract out diffs between this and the primer arrow; avoid repeated code.

use std::ops::Range;

use eframe::{
    egui::{pos2, Align2, Color32, FontFamily, FontId, Pos2, Shape, Stroke, Ui},
    epaint::PathShape,
};

use crate::{
    gui::{
        primer_arrow::{HEIGHT, LABEL_OFFSET, SLANT, STROKE_WIDTH},
        seq_view::NT_WIDTH_PX,
    },
    sequence::{
        Feature, FeatureDirection,
        FeatureDirection::{Forward, Reverse},
        FeatureType,
    },
    util::get_feature_ranges,
    Color,
};
use crate::gui::seq_view::SEQ_ROW_SPACING_PX;

const VERTICAL_OFFSET_FEATURE: f32 = 14.; // Number of pixels above the sequence text.

pub fn draw_features(
    features: &[Feature],
    row_ranges: &[Range<usize>],
    ui: &mut Ui,
    seq_i_to_px_rel: impl Fn(usize) -> Pos2,
) -> Vec<Shape> {
    let mut shapes = Vec::new();

    // todo: Do not run these calcs each time. Cache.
    for feature in features {
        let feature_ranges =
            get_feature_ranges(&(feature.index_range.0..feature.index_range.1), row_ranges);

        let feature_ranges_px: Vec<(Pos2, Pos2)> = feature_ranges
            .iter()
            .map(|r| (seq_i_to_px_rel(r.start), seq_i_to_px_rel(r.end)))
            .collect();

        // todo: PUt back; temp check on compiling.
        shapes.append(&mut feature_seq_overlay(
            &feature_ranges_px,
            feature.feature_type,
            feature.color_override,
            VERTICAL_OFFSET_FEATURE,
            feature.direction,
            &feature.label,
            ui,
        ));
    }
    shapes
}

/// Make a visual indicator on the sequence view for a feature, including primers. For use inside a Frame::canvas.
pub fn feature_seq_overlay(
    feature_ranges_px: &[(Pos2, Pos2)],
    feature_type: FeatureType,
    color_override: Option<Color>,
    vertical_offset: f32,
    direction: FeatureDirection,
    label: &str,
    ui: &mut Ui,
) -> Vec<Shape> {
    if feature_ranges_px.is_empty() {
        return Vec::new();
    }
    let (r, g, b) = match color_override {
        Some(c) => c,
        None => feature_type.color(),
    };
    let stroke = Stroke::new(STROKE_WIDTH, Color32::from_rgb(r, g, b));

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

    for (i, (mut start, mut end)) in feature_ranges_px.iter().enumerate() {
        // Display the overlay centered around the NT letters, vice above, for non-primer features.
        if feature_type != FeatureType::Primer {
            start.y += SEQ_ROW_SPACING_PX / 2. - 2.;
            end.y += SEQ_ROW_SPACING_PX / 2. - 2.;
        }

        let mut top_left = start;
        let mut top_right = pos2(end.x + NT_WIDTH_PX, end.y);
        let bottom_left = pos2(start.x, start.y + HEIGHT);
        let bottom_right = pos2(end.x + NT_WIDTH_PX, end.y + HEIGHT);

        // Add a slant, if applicable.
        match direction {
            Forward => {
                if i + 1 == feature_ranges_px.len() {
                    top_right.x -= SLANT;
                }
            }
            Reverse => {
                if i == 0 {
                    top_left.x += SLANT;
                }
            }
            _ => (),
        }

        result.push(Shape::Path(PathShape::closed_line(
            vec![top_left, bottom_left, bottom_right, top_right],
            stroke,
        )));
    }

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
