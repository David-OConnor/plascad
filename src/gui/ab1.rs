//! Contains code for viewing AB1 sequencing data, e.g. from Sanger sequencing.

use bio_files::SeqRecordAb1;
use copypasta::{ClipboardContext, ClipboardProvider};
use eframe::{
    egui::{
        Align2, Color32, FontFamily, FontId, Frame, Pos2, Rect, RichText, Sense, Shape, Slider,
        Stroke, Ui, pos2, vec2,
    },
    emath::RectTransform,
    epaint::PathShape,
};
use na_seq::{Nucleotide, seq_to_str_lower};

use crate::{
    feature_db_load::find_features,
    file_io::GenericData,
    gui::{BACKGROUND_COLOR, COL_SPACING, ROW_SPACING},
    misc_types::Metadata,
    state::State,
    util::merge_feature_sets,
};

const NT_COLOR: Color32 = Color32::from_rgb(180, 220, 220);
const NT_WIDTH: f32 = 8.; // pixels
const STROKE_WIDTH_PEAK: f32 = 1.;
const PEAK_WIDTH_DIV2: f32 = 2.;

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
    // let width = ui.available_width();
    //
    // let num_nts_disp = width / NT_WIDTH;

    // let scale_factor = width / num_nts as f32;
    i as f32 * NT_WIDTH
}

/// Plot the peaks and confidence values; draw letters
fn plot(data: &SeqRecordAb1, to_screen: &RectTransform, start_i: usize, ui: &mut Ui) -> Vec<Shape> {
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

    let width = ui.available_width();
    let num_nts_disp = width / NT_WIDTH;

    // Display nucleotides and quality values.
    for i_pos in 0..num_nts_disp as usize {
        let i = start_i + i_pos;
        if data.sequence.is_empty() || i > data.sequence.len() - 1 {
            break;
        }

        let nt = data.sequence[i];

        let x_pos = index_to_posit(i_pos, data.sequence.len(), ui);

        if x_pos > ui.available_width() - 15. {
            continue; // todo: QC the details, if you need to_screen here etc.
        }

        let nt_color = nt_color_map(nt);

        result.push(ui.ctx().fonts(|fonts| {
            Shape::text(
                fonts,
                to_screen * pos2(x_pos, nt_y),
                Align2::CENTER_CENTER,
                &nt.to_str_lower(),
                FontId::new(12., FontFamily::Monospace),
                nt_color,
            )
        }));

        // todo: Display quality values below.
        // let quality = data.quality[i]; // todo: Index check.
    }

    // Display data.
    for i_pos in 0..num_nts_disp as usize * 4 {
        let i = start_i * 4 + i_pos;
        if data.data_ch1.is_empty() || i > data.data_ch1.len() - 1 {
            break;
        }

        let x_pos = index_to_posit(i_pos, data.sequence.len(), ui) / 4.;
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
            let top_left = base_pos + vec2(-PEAK_WIDTH_DIV2, ch as f32 * -data_scaler);
            let top_right = base_pos + vec2(PEAK_WIDTH_DIV2, ch as f32 * -data_scaler);
            let bottom_left = base_pos + vec2(-PEAK_WIDTH_DIV2, 0.);
            let bottom_right = base_pos + vec2(PEAK_WIDTH_DIV2, 0.);

            result.push(ui.ctx().fonts(|fonts| {
                // Shape::Path(PathShape::convex_polygon(
                Shape::Path(PathShape::closed_line(
                    vec![
                        to_screen * top_left,
                        to_screen * bottom_left,
                        to_screen * bottom_right,
                        to_screen * top_right,
                    ],
                    // color,
                    stroke,
                ))
            }));
        }
    }

    result
}

pub fn ab1_page(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        let data = &state.ab1_data[state.active];

        ui.heading("AB1 sequencing view");

        ui.add_space(COL_SPACING * 2.);

        if ui.button(RichText::new("🗐 Copy sequence")).on_hover_text("Copy this sequence to the clipboard.").clicked() {
            let mut ctx = ClipboardContext::new().unwrap();
            ctx.set_contents(seq_to_str_lower(&data.sequence)).unwrap();
        }

        ui.add_space(COL_SPACING);

        if ui
            .button("➕ Create data (e.g. Genbank, PCAD etc) as a new tab.")
            .on_hover_text("Create a new non-AB1 data set\
        that may include features, primers, etc. May be edited, and saved to Genbank, PCAD, or SnapGene formats.")
            .clicked()
        {
            // Note: This segment is almost a duplicate of `State::add_tab` and the similar section in `cloning`.
            let plasmid_name = match &state.tabs_open[state.active].path {
                Some(p) => p.file_name().unwrap().to_str().unwrap_or_default(),
                None => "Plasmid from AB1", // This shouldn't happen, I believe.
            }.to_owned().replace(".ab1", "");

            let generic = GenericData {
                seq: data.sequence.clone(),
                metadata: Metadata {
                    plasmid_name,
                    ..Default::default()
                },
                ..Default::default()
            };


            state.generic.push(generic);

            state.portions.push(Default::default());
            state.volatile.push(Default::default());
            state.tabs_open.push(Default::default());
            state.ab1_data.push(Default::default());

            state.active = state.generic.len() - 1;

            // state.sync_seq_related(None);

            // Annotate. Don't add duplicates.
            let features = find_features(&state.get_seq());
            merge_feature_sets(&mut state.generic[state.active].features, &features)

        }
    });
    ui.add_space(ROW_SPACING / 2.);

    let data = &state.ab1_data[state.active];
    let width = ui.available_width();
    let num_nts_disp = width / NT_WIDTH;
    ui.spacing_mut().slider_width = width - 60.;

    let slider_max = {
        let v = data.sequence.len() as isize - num_nts_disp as isize + 1;
        if v > 0 { v as usize } else { 0 }
    };

    ui.add(Slider::new(&mut state.ui.ab1_start_i, 0..=slider_max));

    ui.add_space(ROW_SPACING / 2.);

    let mut shapes = Vec::new();

    Frame::canvas(ui.style())
        .fill(BACKGROUND_COLOR)
        .show(ui, |ui| {
            let data = &state.ab1_data[state.active];

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

            shapes.append(&mut plot(data, &to_screen, state.ui.ab1_start_i, ui));

            ui.painter().extend(shapes);
        });
}
