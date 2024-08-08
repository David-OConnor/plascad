//! This module contains code related to drawing primer arrows in the sequence view.

use std::ops::Range;

use eframe::egui::{Pos2, Shape, Ui};

use crate::{
    gui::feature_overlay,
    primer::{Primer, PrimerData, PrimerDirection},
    sequence::FeatureType,
    util,
};

pub const STROKE_WIDTH: f32 = 2.;

pub const VERTICAL_OFFSET_PRIMER: f32 = 14.; // Number of pixels above the sequence text.
pub const LABEL_OFFSET: f32 = 7.;
pub const HEIGHT: f32 = 16.;
pub const SLANT: f32 = 12.; // slant different, in pixels, for the arrow.

/// Add primer arrows to the display.
pub fn draw_primers(
    primers: &[Primer],
    row_ranges: &[Range<usize>],
    seq_len: usize,
    ui: &mut Ui,
    seq_i_to_px_rel: impl Fn(usize) -> Pos2,
) -> Vec<Shape> {
    let mut shapes = Vec::new();

    for primer in primers {
        let primer_matches = &primer.volatile.matches_seq;

        // todo: Do not run these calcs each time. Cache.
        for (direction, seq_range) in primer_matches {
            // We currently index primers relative to the end they started.
            let seq_range = match direction {
                PrimerDirection::Forward => seq_range.clone(),
                PrimerDirection::Reverse => (seq_len - seq_range.end)..(seq_len - seq_range.start),
            };

            let feature_ranges = util::get_feature_ranges(&seq_range, row_ranges);

            let feature_ranges_px: Vec<(Pos2, Pos2)> = feature_ranges
                .iter()
                .map(|r| (seq_i_to_px_rel(r.start), seq_i_to_px_rel(r.end)))
                .collect();

            let color = match direction {
                PrimerDirection::Forward => (255, 0, 255),
                PrimerDirection::Reverse => (0, 255, 0),
            };

            // todo: PUt back; temp check on compiling.
            shapes.append(&mut feature_overlay::feature_seq_overlay(
                &feature_ranges_px,
                FeatureType::Primer,
                Some(color),
                VERTICAL_OFFSET_PRIMER,
                (*direction).into(),
                &primer.name,
                ui,
            ));
        }
    }
    shapes
}
