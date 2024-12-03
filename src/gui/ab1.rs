//! Contains code for viewing AB1 sequencing data, e.g. from Sanger sequencing.

use eframe::{
    egui::{
        pos2, vec2, Align2, Color32, FontFamily, FontId, Frame, Pos2, Rect, Sense, Shape, Stroke,
        Ui,
    },
    emath::RectTransform,
    epaint::PathShape,
};
use na_seq::Nucleotide;

use crate::{ab1::SeqRecordAb1, gui::BACKGROUND_COLOR};

const NT_COLOR: Color32 = Color32::from_rgb(180, 220, 220);
const NT_WIDTH: f32 = 8.; // pixels
const STROKE_WIDTH_PEAK: f32 = 1.;

// Peak  heights are normallized, so that the maximum value is this.
const PEAK_MAX_HEIGHT: f32 = 120.;

const COLOR_A: Color32 = Color32::from_rgb(20, 220, 20);
const COLOR_C: Color32 = Color32::from_rgb(130, 130, 255);
const COLOR_T: Color32 = Color32::from_rgb(255, 100, 100);
const COLOR_G: Color32 = Color32::from_rgb(200, 200, 200);

/// This mapping is based off conventions in other software.
fn nt_color_map(nt: Nucleotide) -> Color32 {
    match nt {
        Nucleotide::A => COLOR_A,
        Nucleotide::C => COLOR_C,
        Nucleotide::T => COLOR_T,
        Nucleotide::G => COLOR_G,
    }
}

/// Map sequence index to a horizontal pixel.
fn index_to_posit(i: usize, num_nts: usize, ui: &Ui) -> f32 {
    // todo: Cap if width too small for nt len. Varying zoom, etc.
    let width = ui.available_width();

    let num_nts_disp = width / NT_WIDTH;

    // let scale_factor = width / num_nts as f32;
    i as f32 * NT_WIDTH
}

/// Plot the peaks and confidence values; draw letters
fn plot(data: &SeqRecordAb1, to_screen: &RectTransform, ui: &mut Ui) -> Vec<Shape> {
    let mut result = Vec::new();

    let nt_y = 160.;
    let plot_y = nt_y - 20.;

    // todo: Don't run this calc every time.
    let data_scaler = {
        let mut max_peak = 0;
        for i in 0..data.data_ch1.len() {
            if i > 40 {
                // todo: Sloppy performance saver.
                continue;
            }

            let ch1 = data.data_ch1[i];
            let ch2 = data.data_ch2[i];
            let ch3 = data.data_ch3[i];
            let ch4 = data.data_ch4[i];

            for ch in [ch1, ch2, ch3, ch4] {
                if ch > max_peak {
                    max_peak = ch;
                }
            }
        }

        PEAK_MAX_HEIGHT / max_peak as f32
    };

    // Display nucleotides and quality values.
    for (i, nt) in data.sequence.iter().enumerate() {
        if i > 400 {
            // todo: Sloppy performance saver.
            continue;
        }

        let x_pos = index_to_posit(i, data.sequence.len(), ui);

        if x_pos > ui.available_width() - 15. {
            continue; // todo: QC the details, if you need to_screen here etc.
        }

        let nt_color = nt_color_map(*nt);

        result.push(ui.ctx().fonts(|fonts| {
            Shape::text(
                fonts,
                to_screen * pos2(x_pos, nt_y),
                Align2::CENTER_CENTER,
                nt.as_str(),
                FontId::new(12., FontFamily::Monospace),
                nt_color,
            )
        }));

        // todo: Display quality values below.
        // let quality = data.quality[i]; // todo: Index check.
    }

    // Display data.
    for i in 0..data.data_ch1.len() {
        if i > 2_000 {
            // todo: Sloppy performance saver.
            continue;
        }

        let x_pos = index_to_posit(i, data.sequence.len(), ui) / 4.;
        if x_pos > ui.available_width() - 15. {
            continue; // todo: QC the details, if you need to_screen here etc.
        }

        let ch1 = data.data_ch1[i];
        let ch2 = data.data_ch2[i];
        let ch3 = data.data_ch3[i];
        let ch4 = data.data_ch4[i];
        // todo: Index error handling.

        for (ch, color) in [
            (ch1, COLOR_G),
            (ch2, COLOR_A),
            (ch3, COLOR_T),
            (ch4, COLOR_C),
        ] {
            let stroke = Stroke::new(STROKE_WIDTH_PEAK, color);

            // todo: Autoscale  height
            let base_pos = pos2(x_pos, plot_y);

            // todo: These may get too thin.
            let top_left = base_pos + vec2(-NT_WIDTH / 8., ch as f32 * -data_scaler);
            let top_right = base_pos + vec2(NT_WIDTH / 8., ch as f32 * -data_scaler);
            let bottom_left = base_pos + vec2(-NT_WIDTH / 8., 0.);
            let bottom_right = base_pos + vec2(NT_WIDTH / 8., 0.);

            result.push(ui.ctx().fonts(|fonts| {
                Shape::Path(PathShape::convex_polygon(
                    vec![
                        to_screen * top_left,
                        to_screen * bottom_left,
                        to_screen * bottom_right,
                        to_screen * top_right,
                    ],
                    color,
                    stroke,
                ))
            }));
        }
    }

    result
}

pub fn ab1_page(data: &SeqRecordAb1, ui: &mut Ui) {
    ui.heading("AB1 sequencing view");

    let mut shapes = Vec::new();

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

            shapes.append(&mut plot(data, &to_screen, ui));

            ui.painter().extend(shapes);
        });
}
