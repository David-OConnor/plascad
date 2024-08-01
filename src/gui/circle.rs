//! A module for the circular view of a plasmid

use core::f32::consts::TAU;

use eframe::egui::{Align2, Color32, FontFamily, FontId, Frame, Pos2, pos2, Rect, Sense, Shape, Stroke, Ui, vec2};
use eframe::emath::RectTransform;
use eframe::epaint::{CircleShape, PathShape};
use crate::gui::seq_view::{COLOR_SEQ, FONT_SIZE_SEQ, SEQ_ROW_SPACING_PX};
use crate::sequence::{Feature, Nucleotide};
use crate::util::{pixel_to_seq_i, seq_i_to_pixel};

const BACKBONE_COLOR: Color32 = Color32::from_rgb(180, 180, 180);
const BACKBONE_WIDTH: f32 = 10.;

const TICK_COLOR: Color32 = Color32::from_rgb(220, 220, 220);
const TICK_WIDTH: f32 = 4.;

/// Rreturn the angle in radians of a given sequence index.
fn seq_i_to_angle(seq_i: usize, seq_len: usize) -> f32 {
    TAU * seq_i as f32 / seq_len as f32

    // todo: Include mapping to pixels here, or elsewhere?
}

/// Convert an angle, in radians, to pixels. Our origin is at the top, and the direction
/// is clockwise.
fn angle_to_pixel(angle: f32, radius: f32) -> Pos2 {
    // Offset and reverse the angle to reflect an origin of 0.
    let angle = -(angle - TAU/4.);

    pos2(angle.cos() * radius, angle.sin() * radius)
}


pub fn circle_page(seq: &[Nucleotide], features: &[Feature], ui: &mut Ui) {
    let mut shapes = Vec::new();

    Frame::canvas(ui.style())
        .fill(Color32::from_rgb(10, 20, 10))
        .show(ui, |ui| {
            let (mut response, _painter) = {
                // todo: Sort this out to make effective use of the space. Check the examples
                let desired_size = vec2(ui.available_width(), ui.available_height());
                ui.allocate_painter(desired_size, Sense::click())
            };

            // let to_screen =
            //     RectTransform::from_to(Rect::from_x_y_ranges(0.0..=1.0, -1.0..=1.0), rect);

            let to_screen = RectTransform::from_to(
                Rect::from_min_size(Pos2::ZERO, response.rect.size()),
                response.rect,
            );

            let from_screen = to_screen.inverse();

            let center = to_screen * pos2(ui.available_width() / 2., ui.available_height() / 2. + 500.);
            let radius = ui.available_width() * 0.4;

            // Draw the backbone circle
            // center: Pos2, radius: f32, stroke: impl Into<Stroke>
            shapes.push(Shape::Circle(CircleShape::stroke(center, radius, Stroke::new(BACKBONE_WIDTH, BACKBONE_COLOR))));


            let seq_len = seq.len();
 add
            // Draw ticks every 1kbp.
            for i_div_1k in 0..seq_len / 1_000 { // todo: Confirm this always rounds down.
                let i = i_div_1k * 1_000;

                let angle = seq_i_to_angle(i, seq_len);
                let intersection = angle_to_pixel(angle, radius);

                shapes.push(Shape::line_segment(
                    [intersection, pos2(0., 0.)],
                    Stroke::new(TICK_WIDTH, TICK_COLOR),
                ));

            }

            ui.painter().extend(shapes);
        });
}
