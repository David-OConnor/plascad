//! This module is related to drawing features on the sequence view. It is similar to `primer_arrow.rs`.

// todo: Abstract out diffs between this and the primer arrow; avoid repeated code.

use std::ops::{Range, RangeInclusive};

use eframe::{
    egui::{pos2, Align2, Color32, FontFamily, FontId, Pos2, Shape, Stroke, Ui},
    epaint::PathShape,
};

use crate::{
    gui::{
        primer_arrow::{HEIGHT, LABEL_OFFSET, SLANT_DIV2, STROKE_WIDTH},
        seq_view::{SeqViewData, COLOR_CURSOR, NT_WIDTH_PX, SEQ_ROW_SPACING_PX},
    },
    sequence::{
        Feature, FeatureDirection,
        FeatureDirection::{Forward, Reverse},
        FeatureType,
    },
    util::{get_feature_ranges, RangeIncl},
    Color,
};

const VERTICAL_OFFSET_FEATURE: f32 = 18.; // A fudge factor?

/// We include this in this module because visually, it is very similar to the overlay.
/// Note: At least for now, selection uses 1-based indexing.
pub fn draw_selection(selection: RangeIncl, data: &SeqViewData, ui: &mut Ui) -> Vec<Shape> {
    let mut result = Vec::new();

    if selection.start < 1 || selection.end + 1 > data.seq_len {
        eprintln!("Invalid sequence index");
        return result;
    }

    // Todo: Cache this, and only update it if row_ranges change. See what else you can optimize
    // todo in this way.
    let selection_ranges = get_feature_ranges(&selection, &data.row_ranges);

    let selection_ranges_px: Vec<(Pos2, Pos2)> = selection_ranges
        .iter()
        .map(|r| {
            (
                data.seq_i_to_px_rel(*r.start()),
                data.seq_i_to_px_rel(*r.end()),
            )
        })
        .collect();

    result.append(&mut feature_seq_overlay(
        &selection_ranges_px,
        FeatureType::Selection,
        (COLOR_CURSOR.r(), COLOR_CURSOR.g(), COLOR_CURSOR.b()),
        VERTICAL_OFFSET_FEATURE,
        FeatureDirection::None,
        "",
        ui,
    ));

    result
}

pub fn draw_features(features: &[Feature], data: &SeqViewData, ui: &mut Ui) -> Vec<Shape> {
    let mut result = Vec::new();

    for feature in features {
        // Source features generally take up the whole plasmid length.
        // Alternative: Filter by features that take up the whole length.
        if feature.feature_type == FeatureType::Source {
            continue;
        }

        if *feature.range.start() < 1 {
            eprintln!("Invalid sequence index");
            continue; // 0 is invalid, in 1-based indexing, and will underflow.
        }

        // Todo: Cache this, and only update it if row_ranges change. See what else you can optimize
        // todo in this way.
        let feature_ranges = get_feature_ranges(&feature.range, &data.row_ranges);

        let feature_ranges_px: Vec<(Pos2, Pos2)> = feature_ranges
            .iter()
            .map(|r| {
                (
                    data.seq_i_to_px_rel(*r.start()),
                    data.seq_i_to_px_rel(*r.end()),
                )
            })
            .collect();

        let label = if feature.label.is_empty() {
            &feature.feature_type.to_string()
        } else {
            &feature.label
        };

        result.append(&mut feature_seq_overlay(
            &feature_ranges_px,
            feature.feature_type,
            feature.color(),
            VERTICAL_OFFSET_FEATURE,
            feature.direction,
            label,
            ui,
        ));
    }
    result
}

/// Make a visual indicator on the sequence view for a feature, including primers.
/// For use inside a Frame::canvas.
pub fn feature_seq_overlay(
    feature_ranges_px: &[(Pos2, Pos2)],
    feature_type: FeatureType,
    color: Color,
    vertical_offset: f32,
    direction: FeatureDirection,
    label: &str,
    ui: &mut Ui,
) -> Vec<Shape> {
    if feature_ranges_px.is_empty() {
        return Vec::new();
    }
    let (r, g, b) = color;
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

    // Depends on font size.
    let rev_primer_offset = -1.;

    for (i, (mut start, mut end)) in feature_ranges_px.iter().enumerate() {
        // Display the overlay centered around the NT letters, vice above, for non-primer features.
        if feature_type != FeatureType::Primer {
            start.y += SEQ_ROW_SPACING_PX / 2. - 2.;
            end.y += SEQ_ROW_SPACING_PX / 2. - 2.;
        }

        let mut top_left = start;
        let mut top_right = pos2(end.x + NT_WIDTH_PX, end.y);
        let mut bottom_left = pos2(start.x, start.y + HEIGHT);
        let mut bottom_right = pos2(end.x + NT_WIDTH_PX, end.y + HEIGHT);

        // Display reverse primers below the sequence; this vertically mirrors.
        if feature_type == FeatureType::Primer && direction == Reverse {
            top_left.y += 3. * HEIGHT - rev_primer_offset;
            top_right.y += 3. * HEIGHT - rev_primer_offset;
            bottom_left.y += HEIGHT - rev_primer_offset;
            bottom_right.y += HEIGHT - rev_primer_offset;
        }

        // Add a slant, if applicable.
        match direction {
            Forward => {
                if i + 1 == feature_ranges_px.len() {
                    top_right.x -= SLANT_DIV2;
                    bottom_right.x += SLANT_DIV2;
                }
            }
            Reverse => {
                if i == 0 {
                    top_left.x += SLANT_DIV2;
                    bottom_left.x -= SLANT_DIV2;
                }
            }
            _ => (),
        }

        let shape = match feature_type {
            FeatureType::Selection => Shape::Path(PathShape::convex_polygon(
                vec![top_left, bottom_left, bottom_right, top_right],
                stroke.color,
                stroke,
            )),
            _ => Shape::Path(PathShape::closed_line(
                vec![top_left, bottom_left, bottom_right, top_right],
                stroke,
            )),
        };

        result.push(shape);
    }

    // todo: Examine.
    let label_start_x = match direction {
        Forward => feature_ranges_px[0].0.x,
        Reverse => feature_ranges_px[feature_ranges_px.len() - 1].1.x,
        FeatureDirection::None => feature_ranges_px[0].0.x,
    } + LABEL_OFFSET;

    let mut label_pos = match direction {
        Forward => pos2(label_start_x, feature_ranges_px[0].0.y + LABEL_OFFSET),
        Reverse => pos2(
            label_start_x,
            feature_ranges_px[0].0.y + LABEL_OFFSET + v_offset_rev,
        ),
        FeatureDirection::None => pos2(label_start_x, feature_ranges_px[0].0.y + LABEL_OFFSET), // todo: Examine
    };

    let label_align = if feature_type == FeatureType::Primer && direction == Reverse {
        label_pos.y -= 2.;
        Align2::RIGHT_CENTER
    } else {
        Align2::LEFT_CENTER
    };

    let label = ui.ctx().fonts(|fonts| {
        Shape::text(
            fonts,
            label_pos,
            label_align,
            label,
            FontId::new(13., FontFamily::Proportional),
            color_label,
        )
    });

    result.push(label);
    result
}
