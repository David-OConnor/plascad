//! A module for the circular view of a plasmid

use core::f32::consts::TAU;
use std::cmp::min;
use eframe::egui::{Align2, Color32, FontFamily, FontId, Frame, Pos2, pos2, Rect, Sense, Shape, Stroke, Ui, vec2};
use eframe::emath::RectTransform;
use eframe::epaint::{CircleShape, PathShape};
use crate::gui::seq_view::{COLOR_SEQ, FONT_SIZE_SEQ, SEQ_ROW_SPACING_PX};
use crate::sequence::{Feature, Nucleotide};
use crate::util::{pixel_to_seq_i, seq_i_to_pixel};

const BACKBONE_COLOR: Color32 = Color32::from_rgb(180, 180, 180);
const BACKBONE_WIDTH: f32 = 10.;

const TICK_COLOR: Color32 = Color32::from_rgb(180, 220, 220);
const TICK_WIDTH: f32 = 4.;

const TICK_LEN: f32 = 50.; // in pixels.
const TICK_LEN_DIV_2: f32 = TICK_LEN / 2.;
const TICK_LABEL_OFFSET: f32 = 10.;

// Of the available width or height (whichever is lower).
const CIRCLE_SIZE_RATIO: f32 = 0.4;

/// Rreturn the angle in radians of a given sequence index.
fn seq_i_to_angle(seq_i: usize, seq_len: usize) -> f32 {
    TAU * seq_i as f32 / seq_len as f32

    // todo: Include mapping to pixels here, or elsewhere?
}

/// Convert an angle, in radians, to pixels. Our origin is at the top, and the direction
/// is clockwise.
fn angle_to_pixel(angle: f32, radius: f32) -> Pos2 {
    // Offset and reverse the angle to reflect an origin of 0.
    let angle = (angle - TAU/4.);

    pos2(angle.cos() * radius, angle.sin() * radius)
}

/// Draw a tick every 1kbp.
fn draw_ticks(seq_len: usize, center: Pos2, radius: f32, to_screen: &RectTransform,  ui: &mut Ui) -> Vec<Shape> {
    let mut result = Vec::new();

    for i_div_1k in 0..seq_len / 1_000 { // todo: Confirm this always rounds down.
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
            ( point_outer + vec2(-TICK_LABEL_OFFSET, 0.), Align2::RIGHT_CENTER)
        } else {
           ( point_outer + vec2(TICK_LABEL_OFFSET, 0.), Align2::LEFT_CENTER)
        };

        let label = ui.ctx().fonts(|fonts| {
            Shape::text(
                fonts,
                to_screen * label_pt,
                label_align,
                i.to_string(),
                FontId::new(16., FontFamily::Proportional),
                TICK_COLOR
            )
        });

        result.push(label);
    }
    result
}


pub fn circle_page(seq: &[Nucleotide], features: &[Feature], ui: &mut Ui) {
    let mut shapes = Vec::new();

    Frame::canvas(ui.style())
        .fill(Color32::from_rgb(10, 20, 10))
        .show(ui, |ui| {
            let (mut response, _painter) = {
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

            let radius = if width_min { rect_size.x} else { rect_size.y} * CIRCLE_SIZE_RATIO;

            // Draw the backbone circle
            shapes.push(Shape::Circle(CircleShape::stroke(to_screen * center, radius, Stroke::new(BACKBONE_WIDTH, BACKBONE_COLOR))));

            let seq_len = seq.len();

            shapes.append(&mut draw_ticks(seq_len, center, radius, &to_screen, ui));

            ui.painter().extend(shapes);
        });
}
