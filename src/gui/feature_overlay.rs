//! This module is related to drawing features on the sequence view. It is similar to `primer_arrow.rs`.

// todo: Abstract out diffs between this and the primer arrow; avoid repeated code.

use std::ops::Range;

use eframe::{
    egui::{pos2, Align2, Color32, FontFamily, FontId, Pos2, Shape, Stroke, Ui},
    epaint::PathShape,
};

use crate::{
    gui::primer_arrow::primer_arrow,
    sequence::{
        Feature,
        FeatureDirection::{Forward, Reverse},
        FeatureType,
    },
    util::get_feature_ranges,
};

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
        shapes.append(&mut primer_arrow(
            &feature_ranges_px,
            VERTICAL_OFFSET_FEATURE,
            feature.direction,
            &feature.label,
            ui,
        ));
    }
    shapes
}
