//! A module for the circular view of a plasmid

use core::f32::consts::TAU;

use eframe::{
    egui::{
        Align2, Color32, CursorIcon, FontFamily, FontId, Frame, Pos2, Rect, RichText, ScrollArea,
        Sense, Shape, Slider, Stroke, Ui, pos2, vec2,
    },
    emath::RectTransform,
    epaint::{CircleShape, PathShape},
};
use na_seq::{
    restriction_enzyme::{ReMatch, RestrictionEnzyme},
    seq_to_str_lower,
};

use crate::{
    Selection,
    gui::{
        COL_SPACING, COLOR_RE, COLOR_SEQ, PRIMER_FWD_COLOR, ROW_SPACING, SPLIT_SCREEN_MAX_HEIGHT,
        feature_from_index, feature_table::feature_table, get_cursor_text, lin_maps,
        lin_maps::MINI_DISP_NT_LEN, navigation::NAV_BUTTON_COLOR, select_feature,
    },
    misc_types::{Feature, FeatureDirection, FeatureType},
    primer::Primer,
    state::State,
};

const BACKGROUND_COLOR: Color32 = Color32::from_rgb(10, 20, 10);

const BACKBONE_COLOR: Color32 = Color32::from_rgb(180, 180, 180);
const BACKBONE_WIDTH: f32 = 10.;

pub const TICK_COLOR: Color32 = Color32::from_rgb(180, 220, 220);
const TICK_WIDTH: f32 = 2.;
pub const RE_WIDTH: f32 = 2.;
const TICK_SPACING: usize = 500; // Nucleotides between ticks

const TICK_LEN: f32 = 90.; // in pixels.
const TICK_LEN_DIV_2: f32 = TICK_LEN / 2.;
const TICK_LABEL_OFFSET: f32 = 12.;

pub const FEATURE_OUTLINE_COLOR: Color32 = Color32::from_rgb(200, 200, 255);
// const FEATURE_OUTLINE_HIGHLIGHTED: Color32 = Color32::from_rgb(200, 200, 255);
pub const FEATURE_OUTLINE_SELECTED: Color32 = Color32::RED;

const RE_LEN: f32 = 50.; // in pixels.
const RE_LEN_DIV_2: f32 = RE_LEN / 2.;
const RE_LABEL_OFFSET: f32 = 10.;

// We may use per-feature-type widths, but have this for now.
const FEATURE_WIDTH_DEFAULT: f32 = 26.;
pub const FEATURE_STROKE_WIDTH: f32 = 3.;
const PRIMER_WIDTH: f32 = 54.;
pub const PRIMER_STROKE_WIDTH: f32 = 2.;

const TIP_LEN: f32 = 0.03; // Len of arrow tips, in radians
const TIP_WIDTH_RATIO: f32 = 1.5; // Compared to its feature width.

// Radius, comparedd to the available width or height (whichever is lower).
const CIRCLE_SIZE_RATIO: f32 = 0.42;

// The maximum distance the cursor can be from the circle for various tasks like selection, finding
// the cursor position etc.
const SELECTION_MAX_DIST: f32 = 80.;

const CENTER_TEXT_ROW_SPACING: f32 = 20.;

const FEATURE_SLIDER_WIDTH: f32 = 180.;

// We limit each filled concave shape to a circumfrence segment this long, as part of a workaround to EGUI not having a great
// way to draw concave shapes.
const _MAX_ARC_FILL: f32 = 230.;

const VERITICAL_CIRCLE_OFFSET: f32 = 22.; // Useful for leaving room for the zoomed view.

/// These aguments define the circle, and are used in many places in this module.
pub struct CircleData {
    pub seq_len: usize,
    pub center: Pos2,
    pub radius: f32,
    pub to_screen: RectTransform,
    pub from_screen: RectTransform,
    /// This is `from_screen * center`. We store it here to cache.
    pub center_rel: Pos2,
}

impl CircleData {
    /// Automatically populates `center_rel` and `from_screen`.
    fn new(seq_len: usize, center: Pos2, radius: f32, to_screen: RectTransform) -> Self {
        let center_rel = to_screen * center;
        let from_screen = to_screen.inverse();
        Self {
            seq_len,
            center,
            radius,
            to_screen,
            from_screen,
            center_rel,
        }
    }
}

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

/// Draw a tick every 1kbp.
fn draw_ticks(data: &CircleData, ui: &mut Ui) -> Vec<Shape> {
    let mut result = Vec::new();

    for i_div_1k in 0..data.seq_len / TICK_SPACING {
        let i = i_div_1k * TICK_SPACING;

        let angle = seq_i_to_angle(i, data.seq_len);

        let point_inner =
            angle_to_pixel(angle, data.radius - TICK_LEN_DIV_2) + data.center.to_vec2();
        let point_outer =
            angle_to_pixel(angle, data.radius + TICK_LEN_DIV_2) + data.center.to_vec2();

        result.push(Shape::line_segment(
            [data.to_screen * point_inner, data.to_screen * point_outer],
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
                data.to_screen * label_pt,
                label_align,
                i.to_string(),
                FontId::new(16., FontFamily::Proportional),
                TICK_COLOR,
            )
        }));
        result.push(ui.ctx().fonts(|fonts| {
            Shape::text(
                fonts,
                data.to_screen * label_pt,
                label_align,
                i.to_string(),
                FontId::new(16., FontFamily::Proportional),
                TICK_COLOR,
            )
        }));
    }
    result
}

/// Created a filled-in arc. E.g. for features.
fn draw_filled_arc(
    data: &CircleData,
    angle: (f32, f32),
    width: f32,
    fill_color: Color32,
    stroke: Stroke,
) -> Vec<Shape> {
    let mut points_outer = arc_points(data.center_rel, data.radius + width / 2., angle.0, angle.1);
    let mut points_inner = arc_points(data.center_rel, data.radius - width / 2., angle.0, angle.1);

    points_inner.reverse();

    // todo: Confirm this is going clockwise, for nominal performance reasons.
    points_outer.append(&mut points_inner);

    let mut result = Vec::new();
    result.push(Shape::convex_polygon(points_outer, fill_color, stroke));
    // We divide our result into a single line segment, and multiple filled areas. This is to work around
    // EGUI's limitation regarding filling concave shapes.
    // result.push(Shape::closed_line(points_outer, stroke));

    // Draw filled segments. It appears the limitation is based on segment absolute size, vice angular size.
    // No segment will be larger than our threshold.

    // let circum_segment = data.radius * (angle.1 - angle.0);
    // let num_segments = (MAX_ARC_FILL / circum_segment) as usize + 1;
    // let segment_ang_dist = (angle.1 - angle.0) / num_segments as f32;
    //
    // for i in 0..num_segments {
    //     continue;
    //     let ang_seg = (
    //         angle.0 + i as f32 * segment_ang_dist,
    //         angle.0 + (i + 1) as f32 * segment_ang_dist,
    //     );
    //
    //     println!("ANG SEG: {:.4?}. Orig: {:.4?}", ang_seg, angle);
    //     // todo: DRY
    //     let mut points_outer = arc_points(
    //         data.center_rel,
    //         data.radius + width / 2.,
    //         ang_seg.0,
    //         ang_seg.1,
    //     );
    //     let mut points_inner = arc_points(
    //         data.center_rel,
    //         data.radius - width / 2.,
    //         ang_seg.0,
    //         ang_seg.1,
    //     );
    //
    //     points_inner.reverse();
    //
    //     points_outer.append(&mut points_inner);
    //
    //     // result.push(Shape::convex_polygon(points_outer, fill_color, Stroke::NONE));
    //     result.push(Shape::convex_polygon(
    //         points_outer,
    //         Color32::YELLOW,
    //         Stroke::NONE,
    //     ));
    // }

    // Note: We may need to put something like this back, if we use multiple feature widths, feature layers etc.\
    // todo: Perhaps instead of drawing a pie, we draw slimmer slices.
    let points_patch = arc_points(
        data.center_rel,
        data.radius - width / 2. - stroke.width / 2.,
        angle.0,
        angle.1,
    );
    // points_patch.push(data.center);

    result.push(Shape::convex_polygon(
        points_patch,
        BACKGROUND_COLOR,
        Stroke::NONE,
    ));

    result
}

/// This is a fancy way of saying triangle.
fn draw_arrowhead(
    data: &CircleData,
    width: f32,
    angle: (f32, f32),
    color: Color32,
    stroke: Stroke,
) -> Shape {
    let center = data.center.to_vec2();
    let base_outer = data.to_screen * angle_to_pixel(angle.0, data.radius + width / 2.) + center;
    let base_inner = data.to_screen * angle_to_pixel(angle.0, data.radius - width / 2.) + center;
    let tip = data.to_screen * angle_to_pixel(angle.1, data.radius) + center;

    // Points arranged clockwise for performance reasons.
    let points = if angle.1 > angle.0 {
        vec![base_outer, tip, base_inner]
    } else {
        vec![base_inner, tip, base_outer]
    };

    Shape::convex_polygon(points, color, stroke)
}

fn draw_features(
    features: &[Feature],
    data: &CircleData,
    selected: Selection,
    ui: &mut Ui,
) -> Vec<Shape> {
    let mut result = Vec::new();

    for (i, feature) in features.iter().enumerate() {
        // Draw the arc segment.

        // Source features generally take up the whole plasmid length.
        // Alternative: Filter by features that take up the whole length.
        if feature.feature_type == FeatureType::Source {
            continue;
        }
        // todo: Adjust feature, tick etc width (stroke width, and dimensions from cicle) based on window size.

        // todo: Sort out how to handle feature widths. Byy type?
        let feature_width = FEATURE_WIDTH_DEFAULT;
        // todo: Sort out color, as you have a note elsewhere. Type or custom? Type with avail override?

        let (r, g, b) = feature.color();

        let feature_color = Color32::from_rgb(r, g, b);

        let stroke_color = match selected {
            Selection::Feature(j) => {
                if j == i {
                    FEATURE_OUTLINE_SELECTED
                } else {
                    FEATURE_OUTLINE_COLOR
                }
            }
            _ => FEATURE_OUTLINE_COLOR,
        };

        let stroke = Stroke::new(FEATURE_STROKE_WIDTH, stroke_color);

        let angle_start = seq_i_to_angle(feature.range.start, data.seq_len);
        let mut angle_end = seq_i_to_angle(feature.range.end, data.seq_len);

        // This allows for wraps around the origin, without changing other code.
        if angle_end < angle_start {
            angle_end += TAU;
        }

        // We subtract parts from the start or end angle for the arrow tip, if present.
        let angle = match feature.direction {
            FeatureDirection::None => (angle_start, angle_end),
            FeatureDirection::Forward => (angle_start, angle_end - TIP_LEN),
            FeatureDirection::Reverse => (angle_start + TIP_LEN, angle_end),
        };

        result.append(&mut draw_filled_arc(
            data,
            angle,
            feature_width,
            feature_color,
            stroke,
        ));

        // if i == features.len()  - 1 {
        //     // Egui doesn't support concave fills; convex_polygon will spill into the interior concave part.
        //     // Patch this by filling over this with a circle. This is roughly the inner points plus the center point,
        //     // but slightly inwards as not to override the feature edge.
        //
        //     // Draw this after feature bodies, and before arrows, and all other things.
        //
        //     // Note: This single-circle approach will only work for constant  size feature widths.
        //
        //     // This insert location, or something similar, is required.
        //
        //     // todo: WHy is this overriding the arrow heads?
        //     result.push(Shape::Circle(CircleShape::filled(data.center_rel,
        //                                                   data.radius - FEATURE_WIDTH_DEFAULT / 2. - 3. / 2., BACKGROUND_COLOR)
        //
        //     ));
        // }

        // Draw the label.
        let angle_mid = (angle.0 + angle.1) / 2.;

        let point_mid_outer =
            angle_to_pixel(angle_mid, data.radius + feature_width / 2.) + data.center.to_vec2();

        let (mut label_pt, label_align) = if angle_mid > TAU / 2. {
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

        // If towards the very top of bottom, offset the label vertically.
        if angle_mid < TAU / 16. || angle_mid > 15. * TAU / 16. {
            label_pt += vec2(0., -TICK_LABEL_OFFSET);
        } else if angle_mid > 7. * TAU / 16. && angle_mid < 9. * TAU / 16. {
            label_pt += vec2(0., TICK_LABEL_OFFSET);
        }

        // todo: Rotate, and place the label inside the arc segment. Unfortunately, rotating text
        // todo is currently difficult with EGUI.

        // "At the moment this is somewhat tricky; you will need to first produce a ClippedShape with the
        // text you want in it. Then you can use ctx.tesselate() to turn this into a Mesh. Use the
        // rotate method on that mesh. Then you can draw that with the Painter"
        result.push(ui.ctx().fonts(|fonts| {
            Shape::text(
                fonts,
                data.to_screen * label_pt,
                label_align,
                feature.label(),
                FontId::new(14., FontFamily::Proportional),
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
                data,
                feature_width * TIP_WIDTH_RATIO,
                tip_angle,
                feature_color,
                stroke,
            ));
        }
    }

    result
}

/// todo: C+P from draw_features! Build this into the feature one like you did in seq view.
fn draw_primers(
    primers: &[Primer],
    data: &CircleData,
    selected_item: Selection,
    ui: &mut Ui,
) -> Vec<Shape> {
    let mut result = Vec::new();

    let radius_outer = data.radius + PRIMER_WIDTH / 2.;
    let radius_inner = data.radius - PRIMER_WIDTH / 2.;

    for (i, primer) in primers.iter().enumerate() {
        // todo: Do not run these calcs each time. Cache.
        for prim_match in &primer.volatile.matches {
            let angle_start = seq_i_to_angle(prim_match.range.start, data.seq_len);
            let angle_end = seq_i_to_angle(prim_match.range.end, data.seq_len);
            let angle_mid = (angle_start + angle_end) / 2.;

            let point_start_inner =
                angle_to_pixel(angle_start, radius_inner) + data.center.to_vec2();
            let point_start_outer =
                angle_to_pixel(angle_start, radius_outer) + data.center.to_vec2();

            let point_end_inner = angle_to_pixel(angle_end, radius_inner) + data.center.to_vec2();
            let point_end_outer = angle_to_pixel(angle_end, radius_outer) + data.center.to_vec2();

            let point_mid_outer = angle_to_pixel(angle_mid, radius_outer) + data.center.to_vec2();

            let mut outline_color = prim_match.direction.color();

            if let Selection::Primer(sel_i) = selected_item {
                if sel_i == i {
                    outline_color = Color32::RED;
                }
            }

            let stroke = Stroke::new(PRIMER_STROKE_WIDTH, outline_color);

            result.push(Shape::Path(PathShape::line(
                arc_points(data.center_rel, radius_outer, angle_start, angle_end),
                stroke,
            )));
            result.push(Shape::Path(PathShape::line(
                arc_points(data.center_rel, radius_inner, angle_start, angle_end),
                stroke,
            )));

            // Lines connected the inner and outer arcs.
            result.push(Shape::line_segment(
                [
                    data.to_screen * point_start_inner,
                    data.to_screen * point_start_outer,
                ],
                stroke,
            ));
            result.push(Shape::line_segment(
                [
                    data.to_screen * point_end_inner,
                    data.to_screen * point_end_outer,
                ],
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
                    data.to_screen * label_pt,
                    label_align,
                    &primer.name,
                    FontId::new(14., FontFamily::Proportional),
                    stroke.color,
                )
            }));
        }
    }

    result
}

/// Change the selected feature's range from the sequence or map views.
/// todo: Move out of circle.rs, as it's also used in the seq view now.
pub fn feature_range_sliders(state: &mut State, ui: &mut Ui) {
    if let Selection::Feature(feat_i) = &mut state.ui.selected_item {
        ui.spacing_mut().slider_width = FEATURE_SLIDER_WIDTH;
        let seq_len = state.generic[state.active].seq.len();

        if *feat_i >= state.generic[state.active].features.len() {
            eprintln!("Invalid selected feature");
            state.ui.selected_item = Selection::None;
        } else {
            let feature = &mut state.generic[state.active].features[*feat_i];
            // todo: Handle wraps.
            ui.label("Start:");
            ui.add(Slider::new(&mut feature.range.start, 0..=seq_len));

            ui.add_space(COL_SPACING);
            ui.label("End:");
            ui.add(Slider::new(&mut feature.range.end, 0..=seq_len));

            // todo: Don't let end be before start.
            // if feature.range.start > feature.range.end {
            //     feature.range.end = feature.range.start + 1;
            // }
        }
    }
}

fn top_details(state: &mut State, ui: &mut Ui) {
    // todo: A/R
    // display_filters(&mut state.ui, ui);
    // We re-impl display_filters to, to avoid having reading frames. Change A/R.

    let name = if state.ui.re.unique_cutters_only {
        "Single cut sites"
    } else {
        "RE sites"
    }
    .to_owned();

    ui.label(name);
    ui.checkbox(&mut state.ui.seq_visibility.show_res, "");
    ui.add_space(COL_SPACING / 2.);

    ui.label("Features:");
    ui.checkbox(&mut state.ui.seq_visibility.show_features, "");
    ui.add_space(COL_SPACING / 2.);

    ui.label("Primers:");
    ui.checkbox(&mut state.ui.seq_visibility.show_primers, "");
    ui.add_space(COL_SPACING / 2.);

    // Sliders to edit the feature.
    feature_range_sliders(state, ui);

    ui.add_space(COL_SPACING);
    ui.label("Cursor:");
    let cursor_posit_text = get_cursor_text(state.ui.cursor_seq_i, state.get_seq().len());
    ui.heading(cursor_posit_text);
}

/// Find the sequence index under the cursor, if it is over the sequence.
fn find_cursor_i(cursor_pos: Option<(f32, f32)>, data: &CircleData) -> Option<usize> {
    match cursor_pos {
        Some(p) => {
            let pos_rel = data.from_screen * pos2(p.0, p.1);

            let diff = vec2(pos_rel.x - data.center.x, pos_rel.y - data.center.y);
            let cursor_dist = diff.length();

            if (cursor_dist - data.radius).abs() > SELECTION_MAX_DIST {
                return None;
            }

            // Shifted so 0 is at the top.
            let angle = diff.angle() + TAU / 4.;
            Some(angle_to_seq_i(angle, data.seq_len))
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
    data: &CircleData,
    unique_cutters_only: bool,
    sticky_ends_only: bool,
    ui: &mut Ui,
) -> Vec<Shape> {
    let mut result = Vec::new();
    for (i, re_match) in re_matches.iter().enumerate() {
        if re_match.lib_index >= res.len() {
            eprintln!("Invalid RE selected");
            return result;
        }
        let re = &res[re_match.lib_index];

        if (unique_cutters_only && re_match.match_count > 1)
            || (sticky_ends_only && re.makes_blunt_ends())
        {
            continue;
        }

        let cut_i = re_match.seq_index + 1; // to display in the right place.
        let angle = seq_i_to_angle(cut_i + re.cut_after as usize, data.seq_len);

        let point_inner = angle_to_pixel(angle, data.radius - RE_LEN_DIV_2) + data.center.to_vec2();
        let point_outer = angle_to_pixel(angle, data.radius + RE_LEN_DIV_2) + data.center.to_vec2();

        result.push(Shape::line_segment(
            [data.to_screen * point_inner, data.to_screen * point_outer],
            Stroke::new(RE_WIDTH, COLOR_RE),
        ));

        let (mut label_pt, label_align) = if angle > TAU / 2. {
            (
                point_outer + vec2(-RE_LABEL_OFFSET, 0.),
                Align2::RIGHT_CENTER,
            )
        } else {
            (point_outer + vec2(RE_LABEL_OFFSET, 0.), Align2::LEFT_CENTER)
        };

        // Alternate label vertical position, to reduce changes of overlaps.
        if i % 2 == 0 {
            label_pt.y += 22.;
        }

        result.push(ui.ctx().fonts(|fonts| {
            Shape::text(
                fonts,
                data.to_screen * label_pt,
                label_align,
                &re.name,
                FontId::new(16., FontFamily::Proportional),
                COLOR_RE,
            )
        }));
    }
    result
}

/// For drawing feature data in the center of the circle. This may be used for the feature hovered over,
/// or selected.
fn draw_feature_text(feature: &Feature, data: &CircleData, ui: &mut Ui) -> Vec<Shape> {
    let mut result = Vec::new();

    let labels = vec![
        feature.label.clone(),
        feature.location_descrip(data.seq_len),
        feature.feature_type.to_string(),
    ];

    let (r, g, b) = feature.color();
    let color = Color32::from_rgb(r, g, b);

    let mut i = 0; // Rows

    for label in &labels {
        result.push(draw_text(
            label,
            data.to_screen
                * pos2(
                    data.center.x,
                    data.center.y + i as f32 * CENTER_TEXT_ROW_SPACING - 60.,
                ),
            16.,
            color,
            ui,
        )); // slightly below seq name, ui));
        i += 1;
    }

    for note in &feature.notes {
        // Don't let a note overflow. Note: Wrapping would be preferred to this cutoff.
        let max_len = (0.2 * data.radius) as usize; // Note: This depends on font size.
        // Newlines will interfere with our current line-spacing system.
        let text: String = format!("{}: {}", note.0, note.1.replace('\n', " "))
            .chars()
            .take(max_len)
            .collect();

        result.push(ui.ctx().fonts(|fonts| {
            Shape::text(
                fonts,
                data.to_screen
                    * pos2(
                        data.center.x - data.radius * 0.8,
                        data.center.y + i as f32 * CENTER_TEXT_ROW_SPACING - 60.,
                    ),
                Align2::LEFT_CENTER,
                &text,
                FontId::new(13., FontFamily::Proportional),
                TICK_COLOR,
            )
        }));
        i += 1;
    }

    result
}

/// Similar to `draw_feature_text`.
fn draw_primer_text(primer: &Primer, data: &CircleData, ui: &mut Ui) -> Vec<Shape> {
    let mut result = Vec::new();

    let labels = vec![
        primer.name.clone(),
        primer.location_descrip(),
        primer.description.clone().unwrap_or_default(),
        seq_to_str_lower(&primer.sequence),
    ];

    // let color = match primer.direction {
    //     PrimerDirection::Forward => PRIMER_FWD_COLOR,
    //     PrimerDirection::Reverse => PRIMER_REV_COLOR,
    // };
    //
    // todo: Rev or neutral A/R
    let color = PRIMER_FWD_COLOR;

    let mut i = 0; // Rows

    for label in &labels {
        result.push(draw_text(
            label,
            data.to_screen
                * pos2(
                    data.center.x,
                    data.center.y + i as f32 * CENTER_TEXT_ROW_SPACING - 60.,
                ),
            16.,
            color,
            ui,
        )); // slightly below seq name, ui));
        i += 1;
    }

    result
}

/// Draw text in the center of the circle; eg general plasmid information, or information
/// about a feature. This is the selected feature if available; then hovered-over if available;
/// then general plasmid information.
fn draw_center_text(data: &CircleData, state: &mut State, ui: &mut Ui) -> Vec<Shape> {
    let mut result = Vec::new();
    // todo: This nesting is a bit complicated.
    match &state.ui.selected_item {
        Selection::Feature(feat_i) => {
            if state.generic[state.active].features.len() + 1 < *feat_i {
                eprintln!("Invalid selected feature");
            }
            let feature = &state.generic[state.active].features[*feat_i];
            result.append(&mut draw_feature_text(feature, data, ui));
        }
        Selection::Primer(prim_i) => {
            if *prim_i >= state.generic[state.active].primers.len() {
                eprintln!("Invalid primer.");
                return result;
            }
            let primer = &state.generic[state.active].primers[*prim_i];
            result.append(&mut draw_primer_text(primer, data, ui));
        }
        Selection::None => {
            match &state.ui.feature_hover {
                Some(feat_i) => {
                    if *feat_i >= state.generic[state.active].features.len() {
                        eprintln!("Invalid hover feature");
                    }
                    let feature = &state.generic[state.active].features[*feat_i];
                    result.append(&mut draw_feature_text(feature, data, ui));
                }
                None => {
                    // Display a summary of the plasmid
                    result.push(draw_text(
                        &state.generic[state.active].metadata.plasmid_name,
                        data.center_rel,
                        16.,
                        TICK_COLOR,
                        ui,
                    ));
                    result.push(draw_text(
                        &format!("{} bp", data.seq_len),
                        pos2(data.center_rel.x, data.center_rel.y + 20.),
                        13.,
                        TICK_COLOR,
                        ui,
                    ));
                }
            }
        }
    }

    result
}

// todo: Delete
/// Draw a small view of the sequence under the cursor.
fn _draw_mini_seq(data: &CircleData, state: &State, ui: &mut Ui) -> Vec<Shape> {
    let mut result = Vec::new();

    const OFFSET: Pos2 = pos2(4., 6.);

    let seq_full_len = state.get_seq().len();

    if seq_full_len < 10 {
        return result;
    }

    if let Some(i) = state.ui.cursor_seq_i {
        let seq_mini_len = (ui.available_width() / 11.) as usize;

        let mut start = i as isize - seq_mini_len as isize / 2;
        let mut end = i + seq_mini_len / 2;

        if start < 0 {
            start = 0;
        }
        if end + 1 >= seq_full_len && seq_full_len > 0 {
            end = seq_full_len - 1;
        }

        let start = start as usize;

        let seq = &state.get_seq()[start..end];
        let seq_text = seq_to_str_lower(seq);

        result.push(ui.ctx().fonts(|fonts| {
            Shape::text(
                fonts,
                data.to_screen * OFFSET,
                Align2::LEFT_CENTER,
                &seq_text,
                FontId::new(16., FontFamily::Monospace),
                COLOR_SEQ,
            )
        }));
    }

    result
}

pub fn circle_page(state: &mut State, ui: &mut Ui) {
    let mut shapes = Vec::new();

    // todo: ABility to select light mode, and other tools useful for publication.
    if !state.ui.hide_map_feature_editor {
        // Limit the top section height.
        let screen_height = ui.ctx().available_rect().height();
        let half_screen_height = screen_height / SPLIT_SCREEN_MAX_HEIGHT;

        Frame::none().show(ui, |ui| {
            ScrollArea::vertical()
                .max_height(half_screen_height)
                .show(ui, |ui| {
                    feature_table(state, ui);
                });
        });

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

            ui.add_space(COL_SPACING);

            top_details(state, ui);
        });
    }

    ui.ctx()
        .set_cursor_icon(if state.ui.feature_hover.is_some() {
            CursorIcon::PointingHand
        } else {
            CursorIcon::Default
        });

    Frame::canvas(ui.style())
        .fill(BACKGROUND_COLOR)
        .show(ui, |ui| {
            let (response, _painter) = {
                let desired_size = vec2(ui.available_width(), ui.available_height());
                ui.allocate_painter(desired_size, Sense::click())
            };

            let to_screen = RectTransform::from_to(
                // Rect::from_min_size(pos2(0., -VERITICAL_CIRCLE_OFFSET), response.rect.size()),
                Rect::from_min_size(Pos2::ZERO, response.rect.size()),
                response.rect,
            );

            let rect_size = response.rect.size();

            let center = pos2(rect_size.x / 2., rect_size.y / 2. + VERITICAL_CIRCLE_OFFSET);
            let width_min = rect_size.x < rect_size.y;
            let radius = if width_min { rect_size.x } else { rect_size.y } * CIRCLE_SIZE_RATIO;

            let seq_len = state.get_seq().len();

            let data = CircleData::new(seq_len, center, radius, to_screen);

            let prev_cursor_i = state.ui.cursor_seq_i;
            state.ui.cursor_seq_i = find_cursor_i(state.ui.cursor_pos, &data);

            if prev_cursor_i != state.ui.cursor_seq_i {
                state.ui.feature_hover = feature_from_index(
                    &state.ui.cursor_seq_i,
                    &state.generic[state.active].features,
                );
            }

            select_feature(state, &data.from_screen);

            // Draw the backbone circle
            shapes.push(Shape::Circle(CircleShape::stroke(
                data.center_rel,
                radius - BACKBONE_WIDTH / 2.,
                Stroke::new(BACKBONE_WIDTH, BACKBONE_COLOR),
            )));

            // Draw features first, so other items like ticks will be displayed in front of the concave fill circlex.
            if state.ui.seq_visibility.show_features {
                shapes.append(&mut draw_features(
                    &state.generic[state.active].features,
                    &data,
                    state.ui.selected_item,
                    ui,
                ));
            }

            shapes.append(&mut draw_ticks(&data, ui));

            if state.ui.seq_visibility.show_primers {
                shapes.append(&mut draw_primers(
                    &state.generic[state.active].primers,
                    &data,
                    state.ui.selected_item,
                    ui,
                ));
            }

            // tood: Check mark to edit this visibility on the page
            if state.ui.seq_visibility.show_res {
                shapes.append(&mut draw_re_sites(
                    &state.volatile[state.active].restriction_enzyme_matches,
                    &state.restriction_enzyme_lib,
                    &data,
                    state.ui.re.unique_cutters_only,
                    state.ui.re.sticky_ends_only,
                    ui,
                ));
            }

            shapes.append(&mut draw_center_text(&data, state, ui));

            if let Some(cursor_i) = state.ui.cursor_seq_i {
                shapes.append(&mut lin_maps::lin_map_zoomed(
                    &state.generic[state.active],
                    &data.to_screen,
                    cursor_i,
                    MINI_DISP_NT_LEN,
                    state.ui.selected_item,
                    &state.ui,
                    &state.volatile[state.active].restriction_enzyme_matches,
                    &state.restriction_enzyme_lib,
                    ui,
                ));
            }

            ui.painter().extend(shapes);
        });
}
