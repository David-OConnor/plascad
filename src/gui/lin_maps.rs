//! Contains code for our mini view above the circlular map.

use eframe::{
    egui::{
        pos2, vec2, Align2, Color32, FontFamily, FontId, Frame, Pos2, Rect, Sense, Shape, Stroke,
        Ui,
    },
    emath::RectTransform,
};

use crate::{
    gui::{
        circle::{
            CircleData, FEATURE_OUTLINE_COLOR, FEATURE_OUTLINE_SELECTED, FEATURE_STROKE_WIDTH,
            PRIMER_STROKE_WIDTH, RE_WIDTH,
        },
        BACKGROUND_COLOR, COLOR_RE, LINEAR_MAP_HEIGHT,
    },
    primer::{Primer, PrimerDirection},
    restriction_enzyme::{ReMatch, RestrictionEnzyme},
    sequence::Feature,
    util::{map_linear, RangeIncl},
    Selection, State,
};

// How many nucleotides the zoomed-in display at the top of the page represents.
// A smaller value corresponds to a more zoomed-in display.
pub const MINI_DISP_NT_LEN: usize = 400;
// const MINI_DISP_NT_LEN_DIV2: usize = MINI_DISP_NT_LEN / 2;

const CENTER_Y: f32 = 18.;
pub const OFFSET: Pos2 = pos2(4., 6.);
const Y_START: f32 = OFFSET.y + CENTER_Y;

const FEATURE_HEIGHT: f32 = 18.;
const FEATURE_HEIGHT_DIV2: f32 = FEATURE_HEIGHT / 2.;

const PRIMER_HEIGHT: f32 = 26.;
const PRIMER_HEIGHT_DIV2: f32 = PRIMER_HEIGHT / 2.;
const RE_HEIGHT: f32 = 30.;
const RE_HEIGHT_DIV2: f32 = RE_HEIGHT / 2.;

fn draw_features(
    features: &[Feature],
    seq_len: usize,
    to_screen: &RectTransform,
    disp_range: RangeIncl,
    selected_item: Selection,
    index_to_x: impl Fn(usize) -> f32,
    pixel_left: f32,
    pixel_right: f32,
    ui: &mut Ui,
) -> Vec<Shape> {
    let mut result = Vec::new();

    for (i, feature) in features.iter().enumerate() {
        // todo: DRY with teh main view. helper fn.
        let stroke_color = match selected_item {
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

        let mut feature_range = feature.range;

        // Handle wraps around the origin
        if disp_range.end > seq_len {
            feature_range.start += seq_len;
            feature_range.end += seq_len;
        }

        let contains_start = disp_range.contains(feature_range.start);
        let contains_end = disp_range.contains(feature_range.end);

        // todo: Way to not make this a special case?
        let full_size =
            // feature_range.start < disp_range.start && feature_range.end > disp_range.end && disp_range.start < disp_range.end;
            feature_range.start < disp_range.start && feature_range.end > disp_range.end;

        let (r, g, b) = feature.color();
        let color = Color32::from_rgb(r, g, b);

        if contains_start || contains_end || full_size {
            let start_x = if contains_start {
                index_to_x(feature_range.start)
            } else {
                pixel_left
            };

            let end_x = if contains_end {
                index_to_x(feature_range.end)
            } else {
                pixel_right
            };

            // todo: Arrow heads.

            result.push(Shape::convex_polygon(
                vec![
                    to_screen * pos2(start_x, Y_START - FEATURE_HEIGHT_DIV2),
                    to_screen * pos2(end_x, Y_START - FEATURE_HEIGHT_DIV2),
                    to_screen * pos2(end_x, Y_START + FEATURE_HEIGHT_DIV2),
                    to_screen * pos2(start_x, Y_START + FEATURE_HEIGHT_DIV2),
                ],
                color,
                stroke,
            ));
        }

        // Draw the label in the center.  todo: More locations A/R for long features
        let center_i = (feature_range.end + feature_range.start) / 2;

        if disp_range.contains(center_i) {
            let center_x = index_to_x(center_i);

            // Draw the label after the shape.
            let label_pt = pos2(center_x, OFFSET.y + CENTER_Y);

            result.push(ui.ctx().fonts(|fonts| {
                Shape::text(
                    fonts,
                    to_screen * label_pt,
                    Align2::CENTER_CENTER,
                    feature.label(),
                    FontId::new(13., FontFamily::Proportional),
                    Color32::DARK_GREEN,
                )
            }));
        }
    }
    result
}

/// todo: DRY with draw_features (Similar issue to the non-zoomed cirlc.e
fn draw_primers(
    primers: &[Primer],
    seq_len: usize,
    to_screen: &RectTransform,
    disp_range: RangeIncl,
    selected_item: Selection,
    index_to_x: impl Fn(usize) -> f32,
    pixel_left: f32,
    pixel_right: f32,
    ui: &mut Ui,
) -> Vec<Shape> {
    let mut result = Vec::new();

    for (i, primer) in primers.iter().enumerate() {
        for prim_match in &primer.volatile.matches {
            let mut outline_color = prim_match.direction.color();

            if let Selection::Primer(sel_i) = selected_item {
                if sel_i == i {
                    outline_color = Color32::RED;
                }
            }

            let stroke = Stroke::new(PRIMER_STROKE_WIDTH, outline_color);

            let mut prim_range = prim_match.range;

            // Handle wraps around the origin
            if disp_range.end > seq_len {
                prim_range.start += seq_len;
                prim_range.end += seq_len;
            }

            let contains_start = disp_range.contains(prim_range.start);
            let contains_end = disp_range.contains(prim_range.end);

            // todo: Way to not make this a special case?
            let full_size =
                // feature_range.start < disp_range.start && feature_range.end > disp_range.end && disp_range.start < disp_range.end;
                prim_range.start < disp_range.start && prim_range.end > disp_range.end;

            if contains_start || contains_end || full_size {
                let start_x = if contains_start {
                    index_to_x(prim_range.start)
                } else {
                    pixel_left
                };

                let end_x = if contains_end {
                    index_to_x(prim_range.end)
                } else {
                    pixel_right
                };

                // todo: Arrow heads.

                result.push(Shape::convex_polygon(
                    vec![
                        to_screen * pos2(start_x, Y_START - PRIMER_HEIGHT_DIV2),
                        to_screen * pos2(end_x, Y_START - PRIMER_HEIGHT_DIV2),
                        to_screen * pos2(end_x, Y_START + PRIMER_HEIGHT_DIV2),
                        to_screen * pos2(start_x, Y_START + PRIMER_HEIGHT_DIV2),
                    ],
                    Color32::TRANSPARENT,
                    stroke,
                ));
            }

            // Draw the label in the center.  todo: More locations A/R for long features
            let center_i = (prim_range.end + prim_range.start) / 2;

            if disp_range.contains(center_i) {
                let center_x = index_to_x(center_i);

                // Draw the label after the shape.
                let label_v_offset = match prim_match.direction {
                    PrimerDirection::Forward => -18.,
                    PrimerDirection::Reverse => 18.,
                };
                // let label_pt = pos2(center_x, Y_START + PRIMER_HEIGHT_DIV2 + label_v_offset);
                let label_pt = pos2(center_x, OFFSET.y + CENTER_Y + label_v_offset);

                result.push(ui.ctx().fonts(|fonts| {
                    Shape::text(
                        fonts,
                        to_screen * label_pt,
                        Align2::CENTER_CENTER,
                        &primer.name,
                        FontId::new(13., FontFamily::Proportional),
                        outline_color,
                    )
                }));
            }
        }
    }
    result
}

/// Draw RE cut sites through the circle.
/// todo: DRY with tick drawing code.
fn draw_re_sites(
    re_matches: &[ReMatch],
    res: &[RestrictionEnzyme],
    to_screen: &RectTransform,
    index_to_x: impl Fn(usize) -> f32,
    unique_cutters_only: bool,
    ui: &mut Ui,
) -> Vec<Shape> {
    let mut result = Vec::new();
    for (i, re_match) in re_matches.iter().enumerate() {
        if unique_cutters_only && re_match.match_count > 1 {
            continue;
        }

        let cut_i = re_match.seq_index + 1; // to display in the right place.
        let re = &res[re_match.lib_index];

        let point_top = pos2(index_to_x(cut_i), Y_START - RE_HEIGHT_DIV2);
        let point_bottom = pos2(index_to_x(cut_i), Y_START + RE_HEIGHT_DIV2);

        result.push(Shape::line_segment(
            [to_screen * point_bottom, to_screen * point_top],
            Stroke::new(RE_WIDTH, COLOR_RE),
        ));

        let (mut label_pt, label_align) = (point_top + vec2(20., -2.), Align2::LEFT_CENTER);

        // Alternate label vertical position, to reduce changes of overlaps.
        if i % 2 == 0 {
            label_pt.y += RE_HEIGHT + 8.;
        }

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

/// A general purpose linear sequence view, used on several pages.
/// This shows features, primers, and index ticks.
pub fn draw_linear_map(
    state: &State,
    to_screen: &RectTransform,
    index_left: usize,
    index_right: usize,
    show_re_sites: bool,
    active: usize,
    ui: &mut Ui,
) -> Vec<Shape> {
    let mut result = Vec::new();

    let seq_full_len = state.generic[active].seq.len();

    if seq_full_len < 10 {
        return result;
    }

    let pixel_left = OFFSET.x;
    let pixel_right = ui.available_width() - 2. * OFFSET.x;

    let mut disp_range = RangeIncl::new(index_left, index_right);

    // Handle wraps around the origin.
    if disp_range.start > disp_range.end {
        disp_range.end += seq_full_len;
    }

    // Based on the display range, map index to a pixel position.
    let index_to_x = |mut i: usize| {
        // This handles the case when the zoomed-in view is near the top; the left index
        // will be near the end of the sequence, incorrectly calculating the portion-through in
        // the linear map.
        let right = if index_left > index_right {
            if i < index_right {
                // Ie, we are to the right of the origin.
                i += seq_full_len
            }

            index_right + seq_full_len
        } else {
            index_right
        };

        map_linear(
            i as f32,
            (index_left as f32, right as f32),
            (pixel_left, pixel_right),
        )
    };

    result.append(&mut draw_features(
        &state.generic[active].features,
        seq_full_len,
        to_screen,
        disp_range,
        state.ui.selected_item,
        index_to_x,
        pixel_left,
        pixel_right,
        ui,
    ));

    result.append(&mut draw_primers(
        &state.generic[active].primers,
        seq_full_len,
        to_screen,
        disp_range,
        state.ui.selected_item,
        index_to_x,
        pixel_left,
        pixel_right,
        ui,
    ));

    if state.ui.seq_visibility.show_res && show_re_sites {
        result.append(&mut draw_re_sites(
            &state.volatile[active].restriction_enzyme_matches,
            &state.restriction_enzyme_lib,
            to_screen,
            index_to_x,
            state.ui.re.unique_cutters_only,
            ui,
        ));
    }

    result
}

/// Draw a zoomed-in view around the cursor. We use this as the mini view on the circular map page, and on the cloning
/// to show the vector insert point.
/// set `nt_len` to the sequence len if displaying the whole sequence. Set it to a smaller value for a zoomed in view.
pub fn lin_map_zoomed(
    state: &State,
    to_screen: &RectTransform,
    nt_center: usize,
    nt_len: usize,
    active: usize,
    ui: &mut Ui,
) -> Vec<Shape> {
    let mut result = Vec::new();

    let seq_len = state.generic[active].seq.len();

    // Find the bounds of the indices to display; we use them to map indices to pixels for this
    // mini-display.
    // todo: Only if circular.
    let index_left =
        (nt_center as isize - (nt_len / 2) as isize).rem_euclid(seq_len as isize) as usize; // Rust awk % on negative values.
    let index_right = (nt_center + nt_len / 2) % seq_len;

    result.append(&mut draw_linear_map(
        &state,
        to_screen,
        index_left,
        index_right,
        true,
        active,
        ui,
    ));

    result
}

/// Draw a mini sequence display in its own canvas. This displays the entire sequence, and is used on several pages.
/// We use this where we are not drawing on an existing canvas.
pub fn seq_lin_disp(state: &State, ui: &mut Ui, show_re_sites: bool, active: usize) {
    Frame::canvas(ui.style())
        .fill(BACKGROUND_COLOR)
        .show(ui, |ui| {
            let (response, _painter) = {
                let desired_size = vec2(ui.available_width(), LINEAR_MAP_HEIGHT);
                ui.allocate_painter(desired_size, Sense::click())
            };

            let to_screen = RectTransform::from_to(
                Rect::from_min_size(Pos2::ZERO, response.rect.size()),
                response.rect,
            );

            let shapes = draw_linear_map(
                &state,
                &to_screen,
                0,
                state.generic[active].seq.len() - 1,
                show_re_sites,
                active,
                ui,
            );
            ui.painter().extend(shapes);
        });
}
