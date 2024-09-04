//! This module contains code related to drawing primer arrows in the sequence view.

use eframe::egui::{Pos2, Shape, Ui};

use crate::{
    gui::{feature_overlay, seq_view::SeqViewData, PRIMER_FWD_COLOR, PRIMER_REV_COLOR},
    primer::{Primer, PrimerDirection},
    sequence::FeatureType,
    util,
    util::RangeIncl,
    Selection,
};

pub const STROKE_WIDTH: f32 = 2.;

pub const VERTICAL_OFFSET_PRIMER: f32 = 18.; // A fudge factor?
pub const LABEL_OFFSET: f32 = 7.;
pub const HEIGHT: f32 = 16.;
pub const SLANT: f32 = 12.; // slant different, in pixels, for the arrow.
pub const SLANT_DIV2: f32 = SLANT / 2.;

/// Add primer arrows to the display.
pub fn draw_primers(
    primers: &[Primer],
    selected_item: Selection,
    data: &SeqViewData,
    ui: &mut Ui,
) -> Vec<Shape> {
    let mut shapes = Vec::new();

    for (i, primer) in primers.iter().enumerate() {
        let primer_matches = &primer.volatile.matches;

        // todo: Do not run these calcs each time. Cache.
        for prim_match in primer_matches {
            // We currently index primers relative to the end they started.

            // Note: Because if we're displaying above the seq and below, the base of the arrow must match,
            // hence the offset.
            let seq_range = match prim_match.direction {
                PrimerDirection::Forward => {
                    // todo: Getting an underflow, but not sure why yet.
                    let end = if prim_match.range.end > 0 {
                        prim_match.range.end - 1
                    } else {
                        prim_match.range.end
                    };
                    RangeIncl::new(prim_match.range.start, end)
                }

                PrimerDirection::Reverse => {
                    RangeIncl::new(prim_match.range.start + 1, prim_match.range.end)
                }
            };

            let feature_ranges =
                util::get_feature_ranges(&seq_range, &data.row_ranges, data.seq_len);

            let feature_ranges_px: Vec<(Pos2, Pos2)> = feature_ranges
                .iter()
                .map(|r| (data.seq_i_to_px_rel(r.start), data.seq_i_to_px_rel(r.end)))
                .collect();

            let color = match prim_match.direction {
                PrimerDirection::Forward => (
                    PRIMER_FWD_COLOR.r(),
                    PRIMER_FWD_COLOR.g(),
                    PRIMER_FWD_COLOR.b(),
                ),
                PrimerDirection::Reverse => (
                    PRIMER_REV_COLOR.r(),
                    PRIMER_REV_COLOR.g(),
                    PRIMER_REV_COLOR.b(),
                ),
            };

            let selected = match selected_item {
                Selection::Primer(j) => i == j,
                _ => false,
            };

            // todo: PUt back; temp check on compiling.
            shapes.append(&mut feature_overlay::feature_seq_overlay(
                &feature_ranges_px,
                FeatureType::Primer,
                color,
                VERTICAL_OFFSET_PRIMER,
                (prim_match.direction).into(),
                &primer.name,
                selected,
                ui,
            ));
        }
    }
    shapes
}
