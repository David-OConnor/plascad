//! A module for the circular view of a plasmid

use core::f32::consts::TAU;

use eframe::{
    egui::{
        pos2, vec2, Align2, Color32, FontFamily, FontId, Frame, Pos2, Rect, RichText, Sense,
        Shape, Stroke, Ui, Vec2,
    },
    emath::RectTransform,
    epaint::{CircleShape, PathShape},
};

use crate::{
    gui::{
        feature_from_index, features::feature_table, get_cursor_text, navigation::NAV_BUTTON_COLOR,
        primer_arrow::STROKE_WIDTH, seq_view::COLOR_RE, COL_SPACING, ROW_SPACING,
    },
    primer::{Primer, PrimerData, PrimerDirection},
    restriction_enzyme::{ReMatch, RestrictionEnzyme},
    sequence::{Feature, FeatureDirection, FeatureType},
    State,
};

const BACKGROUND_COLOR: Color32 = Color32::from_rgb(10, 20, 10);

const BACKBONE_COLOR: Color32 = Color32::from_rgb(180, 180, 180);
const BACKBONE_WIDTH: f32 = 10.;

const TICK_COLOR: Color32 = Color32::from_rgb(180, 220, 220);
const TICK_WIDTH: f32 = 2.;
const RE_WIDTH: f32 = 2.;
const TICK_SPACING: usize = 500; // Nucleotides between ticks

const FEATURE_OUTLINE_COLOR: Color32 = Color32::from_rgb(200, 200, 255);

const TICK_LEN: f32 = 120.; // in pixels.
const TICK_LEN_DIV_2: f32 = TICK_LEN / 2.;
const TICK_LABEL_OFFSET: f32 = 10.;

const RE_LEN: f32 = 40.; // in pixels.
const RE_LEN_DIV_2: f32 = RE_LEN / 2.;
const RE_LABEL_OFFSET: f32 = 10.;

// We may use per-feature-type widths, but have this for now.
const FEATURE_WIDTH_DEFAULT: f32 = 30.;
const PRIMER_WIDTH: f32 = 54.;
const PRIMER_STROKE_WIDTH: f32 = 2.;

const TIP_LEN: f32 = 0.05; // Len of arrow tips, in radians
const TIP_WIDTH_RATIO: f32 = 1.5; // Compared to its feature width.

// Of the available width or height (whichever is lower).
const CIRCLE_SIZE_RATIO: f32 = 0.4;

// The maximum distance the cursor can be from the circle for various tasks like selection, finding
// the cursor position etc.
const SELECTION_MAX_DIST: f32 = 100.;

// // todo: Experimenting with making a `Mesh`, so we can fill it.
// /// Tessellate a single [`ArcPieShape`] into a [`Mesh`].
// ///
// /// * `arc_pie_shape`: the arc or pie to tessellate.
// /// * `out`: triangles are appended to this.
// pub fn tessellate_arc_pie(shape: &mut Shape, out: &mut Mesh) {
//     // let ArcPieShape {
//     //     center,
//     //     radius,
//     //     start_angle,
//     //     end_angle,
//     //     closed,
//     //     fill,
//     //     stroke,
//     // } = arc_pie_shape;
//
//     if radius <= 0.0
//         || start_angle == end_angle
//         || stroke.width <= 0.0 && (!closed || fill == Color32::TRANSPARENT)
//     {
//         return;
//     }
//
//     if shape.options.coarse_tessellation_culling
//         && !self
//         .clip_rect
//         .expand(radius + stroke.width)
//         .contains(center)
//     {
//         return;
//     }
//
//     // If the arc is a full circle, we can just use the circle function.
//     if (end_angle - start_angle).abs() >= std::f32::consts::TAU {
//         let stroke_color = match stroke.color {
//             ColorMode::Solid(color) => color,
//             ColorMode::UV(callback) => {
//                 // TODO(emilk): Currently, CircleShape does not support PathStroke.
//                 // As a workaround, the stroke color is set to the center color.
//                 // This needs to be revisited once CircleShape gains PathStroke support.
//                 callback(Rect::from_center_size(center, Vec2::splat(radius)), center)
//             }
//         };
//         let stroke = Stroke::new(stroke.width, stroke_color);
//         let circle = CircleShape {
//             center,
//             radius,
//             fill,
//             stroke,
//         };
//         return shape.tessellate_circle(circle, out);
//     }
//
//     shape.scratchpad_path.clear();
//
//     if closed {
//         shape.scratchpad_path
//             .add_pie(center, radius, start_angle, end_angle);
//         shape.scratchpad_path.fill(shape.feathering, fill, out);
//         shape.scratchpad_path
//             .stroke_closed(shape.feathering, &stroke, out);
//     } else {
//         shape.scratchpad_path
//             .add_arc(center, radius, start_angle, end_angle);
//         shape.scratchpad_path
//             .stroke_open(shape.feathering, &stroke, out);
//     }
// }

/// Create points for an arc. Can be used with line_segment to draw the arc.
/// Two of these can be used with convex_polygon to make a filled  arc segment.
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

/// Return the angle in radians of a given sequence index.
fn seq_i_to_angle(seq_i: usize, seq_len: usize) -> f32 {
    TAU * seq_i as f32 / seq_len as f32
    // todo: Include mapping to pixels here, or elsewhere?
}

/// Return the sequence index, corresponding to an angle in radians.
fn angle_to_seq_i(mut angle: f32, seq_len: usize) -> usize {
    if angle < 0. {
        // Our index-finding logic will fail unless the angle is positive.
        angle += TAU;
    }
    (angle * seq_len as f32 / TAU) as usize
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

    for i_div_1k in 0..seq_len / TICK_SPACING {
        let i = i_div_1k * TICK_SPACING;

        let angle = seq_i_to_angle(i, seq_len);

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

/// Created a filled-in arc.
fn draw_filled_arc(
    center: Pos2,
    radius: f32,
    angle: (f32, f32),
    width: f32,
    to_screen: &RectTransform,
    fill_color: Color32,
    stroke: Stroke,
) -> Vec<Shape> {
    let center_screen = to_screen * center;
    let mut points_outer = arc_points(center_screen, radius + width / 2., angle.0, angle.1);
    let mut points_inner = arc_points(center_screen, radius - width / 2., angle.0, angle.1);

    points_inner.reverse();

    // todo: Confirm this is going clockwise, for nominal performance reasons.
    points_outer.append(&mut points_inner);

    let mut result = Vec::new();
    result.push(Shape::convex_polygon(points_outer, fill_color, stroke));

    // Egui doesn't support concave fills; convex_polygon will spill into the interior concave part.
    // Patch this by filling over this with a circle. This is roughly the inner points plus the center point,
    // but slightly inwards as not to override the feature edge.

    let mut points_patch = arc_points(
        center_screen,
        radius - width / 2. - stroke.width / 2.,
        angle.0,
        angle.1,
    );
    points_patch.push(center);

    result.push(Shape::convex_polygon(
        points_patch,
        BACKGROUND_COLOR,
        Stroke::NONE,
    ));

    result
}

/// This is a fancy way of saying triangle.
fn draw_arrowhead(
    center: Vec2,
    radius: f32,
    width: f32,
    angle: (f32, f32),
    color: Color32,
    to_screen: &RectTransform,
) -> Shape {
    let base_outer = to_screen * angle_to_pixel(angle.0, radius + width / 2.) + center;
    let base_inner = to_screen * angle_to_pixel(angle.0, radius - width / 2.) + center;
    let tip = to_screen * angle_to_pixel(angle.1, radius) + center;

    // Points arranged clockwise for performance reasons.
    let points = if angle.1 > angle.0 {
        vec![base_outer, tip, base_inner]
    } else {
        vec![base_inner, tip, base_outer]
    };

    let stroke = Stroke::new(STROKE_WIDTH, FEATURE_OUTLINE_COLOR);

    Shape::convex_polygon(points, color, stroke)
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
        // Draw the arc segment.

        // Source features generally take up the whole plasmid length.
        // Alternative: Filter by features that take up the whole length.
        if feature.feature_type == FeatureType::Source {
            continue;
        }
        // todo: Adjust feature, tick etc width (stroke width, and dimensions from cicle) based on window size.

        // todo: Sort out how to handle feature widths. Byy type?
        let feature_width = FEATURE_WIDTH_DEFAULT;
        let feature_stroke_width = 3.;
        // todo: Sort out color, as you have a note elsewhere. Type or custom? Type with avail override?

        let (r, g, b) = match feature.color_override {
            Some(c) => c,
            None => feature.feature_type.color(),
        };
        let feature_color = Color32::from_rgb(r, g, b);
        let stroke = Stroke::new(feature_stroke_width, FEATURE_OUTLINE_COLOR);

        let angle_start = seq_i_to_angle(feature.index_range.0, seq_len);
        let angle_end = seq_i_to_angle(feature.index_range.1, seq_len);

        // We subtract parts from the start or end angle for the arrow tip, if present.
        let angle = match feature.direction {
            FeatureDirection::None => (angle_start, angle_end),
            FeatureDirection::Forward => (angle_start, angle_end - TIP_LEN),
            FeatureDirection::Reverse => (angle_start + TIP_LEN, angle_end),
        };

        result.append(&mut draw_filled_arc(
            center,
            radius,
            angle,
            feature_width,
            to_screen,
            feature_color,
            stroke,
        ));

        // Draw the label.

        let angle_mid = (angle.0 + angle.1) / 2.;

        let point_mid_outer =
            angle_to_pixel(angle_mid, radius + feature_width / 2.) + center.to_vec2();

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

        // Draw the tip
        if feature.direction != FeatureDirection::None {
            let tip_angle = match feature.direction {
                FeatureDirection::Forward => (angle_end - TIP_LEN, angle_end),
                FeatureDirection::Reverse => (angle_start + TIP_LEN, angle_start),
                _ => unreachable!(),
            };

            result.push(draw_arrowhead(
                center.to_vec2(),
                radius,
                feature_width * TIP_WIDTH_RATIO,
                tip_angle,
                feature_color,
                to_screen,
            ));
        }
    }

    result
}

/// todo: C+P from draw_features! Build this into the feature one like you did in seq view.
fn draw_primers(
    primers: &[Primer],
    seq_len: usize,
    center: Pos2,
    radius: f32,
    to_screen: &RectTransform,
    ui: &mut Ui,
) -> Vec<Shape> {
    let mut result = Vec::new();

    let radius_outer = radius + PRIMER_WIDTH / 2.;
    let radius_inner = radius - PRIMER_WIDTH / 2.;

    for primer in primers {
        let primer_matches = &primer.volatile.matches_seq;

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

            let point_start_inner = angle_to_pixel(angle_start, radius_inner) + center.to_vec2();
            let point_start_outer = angle_to_pixel(angle_start, radius_outer) + center.to_vec2();

            let point_end_inner = angle_to_pixel(angle_end, radius_inner) + center.to_vec2();
            let point_end_outer = angle_to_pixel(angle_end, radius_outer) + center.to_vec2();

            let point_mid_outer = angle_to_pixel(angle_mid, radius_outer) + center.to_vec2();

            // todo: This color code is DRY from primer_arrow. Consolidate.
            let outline_color = match direction {
                PrimerDirection::Forward => Color32::from_rgb(255, 0, 255),
                PrimerDirection::Reverse => Color32::LIGHT_YELLOW,
                // FeatureDirection::None => Color32::GOLD,
            };

            let stroke = Stroke::new(PRIMER_STROKE_WIDTH, outline_color);

            result.push(Shape::Path(PathShape::line(
                arc_points(to_screen * center, radius_outer, angle_start, angle_end),
                stroke,
            )));
            result.push(Shape::Path(PathShape::line(
                arc_points(to_screen * center, radius_inner, angle_start, angle_end),
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
                    &primer.name,
                    FontId::new(16., FontFamily::Proportional),
                    stroke.color,
                )
            }));
        }
    }

    result
}

fn top_details(state: &mut State, ui: &mut Ui) {
    // todo: A/R
    // display_filters(&mut state.ui, ui);
    ui.add_space(COL_SPACING);
    ui.label("Cursor:");
    let cursor_posit_text = get_cursor_text(state.ui.cursor_seq_i, state.generic.seq.len());
    ui.heading(cursor_posit_text);
}

/// Find the sequence index under the cursor, if it is over the sequence.
fn find_cursor_i(
    cursor_pos: Option<(f32, f32)>,
    from_screen: &RectTransform,
    seq_len: usize,
    center: Pos2,
    radius: f32,
) -> Option<usize> {
    match cursor_pos {
        Some(p) => {
            let pos_rel = from_screen * pos2(p.0, p.1);

            let diff = vec2(pos_rel.x - center.x, pos_rel.y - center.y);
            let cursor_dist = diff.length();

            if (cursor_dist - radius).abs() > SELECTION_MAX_DIST {
                return None;
            }

            // Shifted so 0 is at the top.
            let angle = diff.angle() + TAU / 4.;
            Some(angle_to_seq_i(angle, seq_len))
        }
        None => None,
    }
}

/// Helper fn.
fn draw_text(text: &str, pos: Pos2, font_size: f32, color: Color32, ui: &mut Ui) -> Shape {
    ui.ctx().fonts(|fonts| {
        Shape::text(
            fonts,
            pos,
            Align2::CENTER_CENTER,
            text,
            FontId::new(font_size, FontFamily::Proportional),
            color,
        )
    })
}

/// Draw RE cut sites through the circle.
/// todo: DRY with tick drawing code.
fn draw_re_sites(
    re_matches: &[ReMatch],
    res: &[RestrictionEnzyme],
    seq_len: usize,
    center: Pos2,
    radius: f32,
    to_screen: &RectTransform,
    ui: &mut Ui,
) -> Vec<Shape> {
    let mut result = Vec::new();
    for re_match in re_matches {
        let cut_i = re_match.seq_index + 1; // to display in the right place.
        let re = &res[re_match.lib_index];
        let angle = seq_i_to_angle(cut_i + re.cut_after as usize, seq_len);

        let point_inner = angle_to_pixel(angle, radius - RE_LEN_DIV_2) + center.to_vec2();
        let point_outer = angle_to_pixel(angle, radius + RE_LEN_DIV_2) + center.to_vec2();

        result.push(Shape::line_segment(
            [to_screen * point_inner, to_screen * point_outer],
            Stroke::new(RE_WIDTH, COLOR_RE),
        ));

        let (label_pt, label_align) = if angle > TAU / 2. {
            (
                point_outer + vec2(-RE_LABEL_OFFSET, 0.),
                Align2::RIGHT_CENTER,
            )
        } else {
            (
                point_outer + vec2(RE_LABEL_OFFSET, 0.),
                Align2::LEFT_CENTER,
            )
        };

        result.push(ui.ctx().fonts(|fonts| {
            Shape::text(
                fonts,
                to_screen * label_pt,
                label_align,
               &re.name,
                FontId::new(16., FontFamily::Proportional),
                COLOR_RE,
            )
        }));
    }
    result
}

/// Draw text in the center of the circle; eg general plasmid information, or information
/// about a feature.
fn draw_center_text(
    center: Pos2,
    to_screen: &RectTransform,
    seq_len: usize,
    state: &mut State,
    ui: &mut Ui,
) -> Vec<Shape> {
    let mut result = Vec::new();
    // todo: Separate function for center label too if it becomes too complicatged.

    match &state.ui.feature_hover {
        Some(i) => {
            if state.generic.features.len() + 1 < *i {
                eprintln!("Invalid hover feature");
            }
            let feature = &state.generic.features[*i];

            let labels = vec![
                feature.label.clone(),
                format!("{}..{}", feature.index_range.0, feature.index_range.1),
                feature.feature_type.to_string(),
            ];

            let row_spacing = 20.;
            let mut i = 0; // Rows
            for label in &labels {
                result.push(draw_text(
                    label,
                    to_screen * pos2(center.x, center.y + i as f32 * row_spacing - 60.),
                    16.,
                    Color32::WHITE,
                    ui,
                )); // slightly below seq name, ui));
                i += 1;
            }

            i += 1;

            for note in &feature.notes {
                result.push(draw_text(
                    &format!("{}: {}", note.0, note.1),
                    to_screen * pos2(center.x, center.y + i as f32 * row_spacing - 60.),
                    13.,
                    TICK_COLOR,
                    ui,
                ));
                i += 1;
            }

            // todo: COlor-code etc?
        }
        None => {
            // Display a summary of the plasmid
            result.push(draw_text(
                &state.generic.metadata.plasmid_name,
                to_screen * center,
                16.,
                TICK_COLOR,
                ui,
            ));
            result.push(draw_text(
                &format!("{seq_len} bp"),
                to_screen * pos2(center.x, center.y + 20.),
                13.,
                TICK_COLOR,
                ui,
            ));
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
            top_details(state, ui);
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
            top_details(state, ui);
        });
    }

    Frame::canvas(ui.style())
        .fill(BACKGROUND_COLOR)
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

            // todo: Cache to_screen * center etc?

            let prev_cursor_i = state.ui.cursor_seq_i;
            state.ui.cursor_seq_i = find_cursor_i(
                state.ui.cursor_pos,
                &from_screen,
                state.generic.seq.len(),
                center,
                radius,
            );

            if prev_cursor_i != state.ui.cursor_seq_i {
                state.ui.feature_hover = None;
                // todo: Consider cacheing this, instead of running each renderx.
                // todo: You may not need the state.ui hover_feature i: You can probably use a local ref here.
                state.ui.feature_hover =
                    feature_from_index(&state.ui.cursor_seq_i, &state.generic.features);
            }

            // Draw the backbone circle
            shapes.push(Shape::Circle(CircleShape::stroke(
                to_screen * center,
                radius,
                Stroke::new(BACKBONE_WIDTH, BACKBONE_COLOR),
            )));

            let seq_len = state.generic.seq.len();

            shapes.append(&mut draw_ticks(seq_len, center, radius, &to_screen, ui));
            shapes.append(&mut draw_features(
                &state.generic.features,
                seq_len,
                center,
                radius,
                &to_screen,
                ui,
            ));

            shapes.append(&mut draw_primers(
                &state.generic.primers,
                seq_len,
                center,
                radius,
                &to_screen,
                ui,
            ));

            shapes.append(&mut draw_re_sites(
                &state.volatile.restriction_enzyme_sites,
                &state.restriction_enzyme_lib,
                seq_len,
                center,
                radius,
                &to_screen,
                ui,
            ));

            shapes.append(&mut draw_center_text(
                center, &to_screen, seq_len, state, ui,
            ));

            ui.painter().extend(shapes);
        });
}
