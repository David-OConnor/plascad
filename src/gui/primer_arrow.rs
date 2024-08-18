//! This module contains code related to drawing primer arrows in the sequence view.

use eframe::egui::{Pos2, Shape, Ui};

use crate::{
    gui::{feature_overlay, seq_view::SeqViewData},
    primer::{Primer, PrimerDirection},
    sequence::FeatureType,
    util,
    util::RangeIncl,
};

pub const STROKE_WIDTH: f32 = 2.;

pub const VERTICAL_OFFSET_PRIMER: f32 = 18.; // A fudge factor?
pub const LABEL_OFFSET: f32 = 7.;
pub const HEIGHT: f32 = 16.;
pub const SLANT: f32 = 12.; // slant different, in pixels, for the arrow.
pub const SLANT_DIV2: f32 = SLANT / 2.;

/// Add primer arrows to the display.
pub fn draw_primers(primers: &[Primer], data: &SeqViewData, ui: &mut Ui) -> Vec<Shape> {
    let mut shapes = Vec::new();

    for primer in primers {
        let primer_matches = &primer.volatile.matches_seq;

        // todo: Do not run these calcs each time. Cache.
        for (direction, seq_range) in primer_matches {
            // We currently index primers relative to the end they started.
            let seq_range = match direction {
                PrimerDirection::Forward => seq_range.clone(),
                PrimerDirection::Reverse => {
                    RangeIncl::new(data.seq_len - seq_range.end, data.seq_len - seq_range.start)
                }
            };

            let feature_ranges = util::get_feature_ranges(&seq_range, &data.row_ranges, data.seq_len);

            let feature_ranges_px: Vec<(Pos2, Pos2)> = feature_ranges
                .iter()
                .map(|r| (data.seq_i_to_px_rel(r.start), data.seq_i_to_px_rel(r.end)))
                .collect();

            let color = match direction {
                PrimerDirection::Forward => (255, 0, 255),
                PrimerDirection::Reverse => (0, 255, 0),
            };

            // todo: PUt back; temp check on compiling.
            shapes.append(&mut feature_overlay::feature_seq_overlay(
                &feature_ranges_px,
                FeatureType::Primer,
                color,
                VERTICAL_OFFSET_PRIMER,
                (*direction).into(),
                &primer.name,
                false,
                ui,
            ));
        }
    }
    shapes
}
