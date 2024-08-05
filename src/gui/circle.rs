//! A module for the circular view of a plasmid

use core::f32::consts::TAU;
use std::cmp::min;

use eframe::{
    egui::{
        pos2, vec2, Align2, Color32, FontFamily, FontId, Frame, Pos2, Rect, RichText, Sense, Shape,
        Stroke, Ui,
    },
    emath::RectTransform,
    epaint::{CircleShape, PathShape},
};

use crate::{
    gui::{
        features::feature_table,
        navigation::NAV_BUTTON_COLOR,
        seq_view::{display_filters, COLOR_SEQ, FONT_SIZE_SEQ, SEQ_ROW_SPACING_PX},
        COL_SPACING, ROW_SPACING,
    },
    primer::{PrimerData, PrimerDirection},
    sequence::{Feature, Nucleotide},
    util::{pixel_to_seq_i, seq_i_to_pixel},
    State,
};

const BACKBONE_COLOR: Color32 = Color32::from_rgb(180, 180, 180);
const BACKBONE_WIDTH: f32 = 10.;

const TICK_COLOR: Color32 = Color32::from_rgb(180, 220, 220);
const TICK_WIDTH: f32 = 3.;

const TICK_LEN: f32 = 46.; // in pixels.
const TICK_LEN_DIV_2: f32 = TICK_LEN / 2.;
const TICK_LABEL_OFFSET: f32 = 10.;

// Of the available width or height (whichever is lower).
const CIRCLE_SIZE_RATIO: f32 = 0.4;

/// Create points for an arc. Can be used with line_segment to draw the arc.
// Adapted from PR https://github.com/emilk/egui/pull/4836/files
fn arc_points(center: Pos2, radius: f32, start_angle: f32, end_angle: f32) -> Vec<Pos2> {
    let num_segs = if radius <= 2.0 {
        8
    } else if radius <= 5.0 {
        16
    } else if radius < 18.0 {
        32
    } else if radius < 50.0 {
        64
    } else {
        128
    };

    let angle = (end_angle - start_angle).clamp(-TAU + f32::EPSILON, TAU - f32::EPSILON);
    let mut points = Vec::with_capacity(num_segs + 3);
    let step = angle / num_segs as f32;

    for i in 0..=num_segs {
        let a = start_angle + step * i as f32;
        points.push(angle_to_pixel(a, radius) + center.to_vec2());
    }

    points
}

/// Rreturn the angle in radians of a given sequence index.
fn seq_i_to_angle(seq_i: usize, seq_len: usize) -> f32 {
    TAU * seq_i as f32 / seq_len as f32

    // todo: Include mapping to pixels here, or elsewhere?
}

/// Convert an angle, in radians, to pixels. Our origin is at the top, and the direction
/// is clockwise.
fn angle_to_pixel(angle: f32, radius: f32) -> Pos2 {
    // Offset and reverse the angle to reflect an origin of 0.
    let angle = angle - TAU / 4.;

    pos2(angle.cos() * radius, angle.sin() * radius)
}

// /// Convenience struct that groups related parameters required for drawing most circle items.
// struct CircleConfig {
//     seq_len: usize,
//     center: Pos2,
//     radius: f32,
//     to_screen: RectTransform,
// }

/// Draw a tick every 1kbp.
fn draw_ticks(
    seq_len: usize,
    center: Pos2,
    radius: f32,
    to_screen: &RectTransform,
    ui: &mut Ui,
) -> Vec<Shape> {
    let mut result = Vec::new();

    for i_div_1k in 0..seq_len / 1_000 {
        // todo: Confirm this always rounds down.
        let i = i_div_1k * 1_000;

        let angle = seq_i_to_angle(i, seq_len);
        // let intersection = angle_to_pixel(angle, radius) + center.to_vec2();

        let point_inner = angle_to_pixel(angle, radius - TICK_LEN_DIV_2) + center.to_vec2();
        let point_outer = angle_to_pixel(angle, radius + TICK_LEN_DIV_2) + center.to_vec2();

        result.push(Shape::line_segment(
            [to_screen * point_inner, to_screen * point_outer],
            Stroke::new(TICK_WIDTH, TICK_COLOR),
        ));

        let (label_pt, label_align) = if angle > TAU / 2. {
            (
                point_outer + vec2(-TICK_LABEL_OFFSET, 0.),
                Align2::RIGHT_CENTER,
            )
        } else {
            (
                point_outer + vec2(TICK_LABEL_OFFSET, 0.),
                Align2::LEFT_CENTER,
            )
        };

        result.push(ui.ctx().fonts(|fonts| {
            Shape::text(
                fonts,
                to_screen * label_pt,
                label_align,
                i.to_string(),
                FontId::new(16., FontFamily::Proportional),
                TICK_COLOR,
            )
        }));
    }
    result
}

fn draw_features(
    features: &[Feature],
    seq_len: usize,
    center: Pos2,
    radius: f32,
    to_screen: &RectTransform,
    ui: &mut Ui,
) -> Vec<Shape> {
    let mut result = Vec::new();

    for feature in features {
        let angle_start = seq_i_to_angle(feature.index_range.0, seq_len);
        let angle_end = seq_i_to_angle(feature.index_range.1, seq_len);
        let angle_mid = (angle_start + angle_end) / 2.;

        // todo: Adjust feature, tick etc width (stroke width, and dimensions from cicle) based on window size.

        // todo: Sort out how to handle feature widths. Byy type?
        let feature_width = 40.;
        let feature_stroke_width = 3.;
        // todo: Sort out color, as you have a note elsewhere. Type or custom? Type with avail override?

        let point_start_inner =
            angle_to_pixel(angle_start, radius - feature_width / 2.) + center.to_vec2();
        let point_start_outer =
            angle_to_pixel(angle_start, radius + feature_width / 2.) + center.to_vec2();

        let point_end_inner =
            angle_to_pixel(angle_end, radius - feature_width / 2.) + center.to_vec2();
        let point_end_outer =
            angle_to_pixel(angle_end, radius + feature_width / 2.) + center.to_vec2();

        let point_mid_outer =
            angle_to_pixel(angle_mid, radius + feature_width / 2.) + center.to_vec2();

        let (r, g, b) = match feature.color_override {
            Some(c) => c,
            None => feature.feature_type.color(),
        };
        let stroke = Stroke::new(feature_stroke_width, Color32::from_rgb(r, g, b));

        result.push(Shape::Path(PathShape::line(
            arc_points(
                to_screen * center,
                radius + feature_width / 2.,
                angle_start,
                angle_end,
            ),
            stroke,
        )));
        result.push(Shape::Path(PathShape::line(
            arc_points(
                to_screen * center,
                radius - feature_width / 2.,
                angle_start,
                angle_end,
            ),
            stroke,
        )));

        // Lines connected the inner and outer arcs.
        result.push(Shape::line_segment(
            [to_screen * point_start_inner, to_screen * point_start_outer],
            stroke,
        ));
        result.push(Shape::line_segment(
            [to_screen * point_end_inner, to_screen * point_end_outer],
            stroke,
        ));

        // todo: A/R

        let (label_pt, label_align) = if angle_mid > TAU / 2. {
            (
                point_mid_outer + vec2(-TICK_LABEL_OFFSET, 0.),
                Align2::RIGHT_CENTER,
            )
        } else {
            (
                point_mid_outer + vec2(TICK_LABEL_OFFSET, 0.),
                Align2::LEFT_CENTER,
            )
        };

        result.push(ui.ctx().fonts(|fonts| {
            Shape::text(
                fonts,
                to_screen * label_pt,
                label_align,
                &feature.label,
                FontId::new(16., FontFamily::Proportional),
                stroke.color,
            )
        }));
    }

    result
}

/// todo: C+P from draw_features! Build this into the feature one like you did in seq view.
fn draw_primers(
    primer_data: &[PrimerData],
    seq_len: usize,
    center: Pos2,
    radius: f32,
    to_screen: &RectTransform,
    ui: &mut Ui,
) -> Vec<Shape> {
    let mut result = Vec::new();

    for prim_data in primer_data {
        let primer_matches = &prim_data.matches_seq;

        // todo: Do not run these calcs each time. Cache.
        for (direction, seq_range) in primer_matches {
            // We currently index primers relative to the end they started.
            let seq_range = match direction {
                PrimerDirection::Forward => seq_range.clone(),
                PrimerDirection::Reverse => (seq_len - seq_range.end)..(seq_len - seq_range.start),
            };

            let angle_start = seq_i_to_angle(seq_range.start, seq_len);
            let angle_end = seq_i_to_angle(seq_range.end, seq_len);
            let angle_mid = (angle_start + angle_end) / 2.;

            // todo: Adjust feature, tick etc width (stroke width, and dimensions from cicle) based on window size.

            // todo: Sort out how to handle feature widths. Byy type?
            let feature_width = 40.;
            let feature_stroke_width = 3.;
            // todo: Sort out color, as you have a note elsewhere. Type or custom? Type with avail override?

            let point_start_inner =
                angle_to_pixel(angle_start, radius - feature_width / 2.) + center.to_vec2();
            let point_start_outer =
                angle_to_pixel(angle_start, radius + feature_width / 2.) + center.to_vec2();

            let point_end_inner =
                angle_to_pixel(angle_end, radius - feature_width / 2.) + center.to_vec2();
            let point_end_outer =
                angle_to_pixel(angle_end, radius + feature_width / 2.) + center.to_vec2();

            let point_mid_outer =
                angle_to_pixel(angle_mid, radius + feature_width / 2.) + center.to_vec2();

            // todo: This color code is DRY from primer_arrow. Consolidate.
            let outline_color = match direction {
                PrimerDirection::Forward => Color32::from_rgb(255, 0, 255),
                PrimerDirection::Reverse => Color32::LIGHT_YELLOW,
                // FeatureDirection::None => Color32::GOLD,
            };

            let stroke = Stroke::new(feature_stroke_width, outline_color);

            result.push(Shape::Path(PathShape::line(
                arc_points(
                    to_screen * center,
                    radius + feature_width / 2.,
                    angle_start,
                    angle_end,
                ),
                stroke,
            )));
            result.push(Shape::Path(PathShape::line(
                arc_points(
                    to_screen * center,
                    radius - feature_width / 2.,
                    angle_start,
                    angle_end,
                ),
                stroke,
            )));

            // Lines connected the inner and outer arcs.
            result.push(Shape::line_segment(
                [to_screen * point_start_inner, to_screen * point_start_outer],
                stroke,
            ));
            result.push(Shape::line_segment(
                [to_screen * point_end_inner, to_screen * point_end_outer],
                stroke,
            ));

            // todo: A/R

            let (label_pt, label_align) = if angle_mid > TAU / 2. {
                (
                    point_mid_outer + vec2(-TICK_LABEL_OFFSET, 0.),
                    Align2::RIGHT_CENTER,
                )
            } else {
                (
                    point_mid_outer + vec2(TICK_LABEL_OFFSET, 0.),
                    Align2::LEFT_CENTER,
                )
            };

            result.push(ui.ctx().fonts(|fonts| {
                Shape::text(
                    fonts,
                    to_screen * label_pt,
                    label_align,
                    &prim_data.primer.description,
                    FontId::new(16., FontFamily::Proportional),
                    stroke.color,
                )
            }));
        }
    }

    result
}

pub fn circle_page(state: &mut State, ui: &mut Ui) {
    let mut shapes = Vec::new();

    // todo: ABility to select light mode, and other tools useful for publication.

    if !state.ui.hide_map_feature_editor {
        feature_table(state, ui);

        ui.add_space(ROW_SPACING / 2.);

        ui.horizontal(|ui| {
            if ui
                .button(RichText::new("Hide editor").background_color(NAV_BUTTON_COLOR))
                .clicked()
            {
                state.ui.hide_map_feature_editor = true;
            }

            // todo: A/R
            // ui.add_space(COL_SPACING);
            // display_filters(&mut state.ui, ui);
        });

        ui.add_space(ROW_SPACING / 2.);
    } else {
        ui.horizontal(|ui| {
            if ui
                .button(RichText::new("Show feature editor").background_color(NAV_BUTTON_COLOR))
                .clicked()
            {
                state.ui.hide_map_feature_editor = false;
            }

            // todo: A/R
            // ui.add_space(COL_SPACING);
            // display_filters(&mut state.ui, ui);
        });
    }

    Frame::canvas(ui.style())
        .fill(Color32::from_rgb(10, 20, 10))
        .show(ui, |ui| {
            let (response, _painter) = {
                // todo: Sort this out to make effective use of the space. Check the examples

                // todo: avail height showing 0.
                let desired_size = vec2(ui.available_width(), ui.available_height());

                // let (_id, rect) = ui.allocate_space(desired_size);

                ui.allocate_painter(desired_size, Sense::click())
            };

            let to_screen = RectTransform::from_to(
                Rect::from_min_size(Pos2::ZERO, response.rect.size()),
                response.rect,
            );

            let from_screen = to_screen.inverse();

            let rect_size = response.rect.size();

            let center = pos2(rect_size.x / 2., rect_size.y / 2.);

            let width_min = rect_size.x < rect_size.y;

            let radius = if width_min { rect_size.x } else { rect_size.y } * CIRCLE_SIZE_RATIO;

            // Draw the backbone circle
            shapes.push(Shape::Circle(CircleShape::stroke(
                to_screen * center,
                radius,
                Stroke::new(BACKBONE_WIDTH, BACKBONE_COLOR),
            )));

            let seq_len = state.seq.len();

            shapes.append(&mut draw_ticks(seq_len, center, radius, &to_screen, ui));
            shapes.append(&mut draw_features(
                &state.features,
                seq_len,
                center,
                radius,
                &to_screen,
                ui,
            ));

            shapes.append(&mut draw_primers(
                &state.primer_data,
                seq_len,
                center,
                radius,
                &to_screen,
                ui,
            ));

            // todo: Separate function for center label too if it becomes too complicatged.
            shapes.push(ui.ctx().fonts(|fonts| {
                Shape::text(
                    fonts,
                    to_screen * center,
                    Align2::CENTER_CENTER,
                    &state.plasmid_name,
                    FontId::new(16., FontFamily::Proportional),
                    TICK_COLOR,
                )
            }));

            shapes.push(ui.ctx().fonts(|fonts| {
                Shape::text(
                    fonts,
                    to_screen * pos2(center.x, center.y + 20.), // slightly below seq name
                    Align2::CENTER_CENTER,
                    format!("{seq_len}bp"),
                    FontId::new(13., FontFamily::Proportional),
                    TICK_COLOR,
                )
            }));

            ui.painter().extend(shapes);
        });
}
