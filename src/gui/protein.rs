use eframe::{
    egui::{
        pos2, vec2, Align2, Color32, FontFamily, FontId, Frame, Pos2, Rect, RichText, Sense, Shape,
        Stroke, Ui,
    },
    emath::RectTransform,
    epaint::PathShape,
};

use crate::{
    amino_acids::{AaIdent, AminoAcid},
    gui::{circle::TICK_COLOR, seq_view::BACKGROUND_COLOR, COL_SPACING, ROW_SPACING},
    sequence::FeatureType,
    State,
};

const COLOR_PROT_SEQ: Color32 = Color32::from_rgb(255, 100, 200);
const COLOR_PRE_POST_CODING_SEQ: Color32 = Color32::from_rgb(100, 255, 200);
const FONT_SIZE_SEQ: f32 = 14.;

const CHART_HEIGHT: f32 = 200.;
const CHART_LINE_WIDTH: f32 = 2.;
const CHART_LINE_COLOR: Color32 = Color32::from_rgb(255, 100, 100);

// todo: Color-code AAs, start/stop codons etc.

// todo: Eval how cacheing and state is handled.

/// Convert an AA sequence to an ident string.
fn make_aa_text(seq: &[AminoAcid], aa_ident_disp: AaIdent) -> String {
    let mut result = String::new();
    for aa in seq {
        let aa_str = match aa_ident_disp {
            AaIdent::ThreeLetters => &format!("{}  ", aa.ident_3_letter()),
            AaIdent::OneLetter => &format!("{}  ", aa.ident_single_letter()),
        };
        result.push_str(aa_str);
    }

    result
}

/// Plot a hydrophobicity line plot.
fn hydrophobicity_chart(data: &Vec<(usize, f32)>, ui: &mut Ui) {
    let stroke = Stroke::new(CHART_LINE_WIDTH, CHART_LINE_COLOR);

    Frame::canvas(ui.style())
        .fill(BACKGROUND_COLOR)
        .show(ui, |ui| {
            let width = ui.available_width();
            let (response, _painter) = {
                let desired_size = vec2(width, CHART_HEIGHT);
                // ui.allocate_painter(desired_size, Sense::click())
                ui.allocate_painter(desired_size, Sense::click_and_drag())
            };

            let to_screen = RectTransform::from_to(
                Rect::from_min_size(Pos2::ZERO, response.rect.size()),
                response.rect,
            );
            // let from_screen = to_screen.inverse();

            const MAX_VAL: f32 = 6.;

            let mut points = Vec::new();
            let num_pts = data.len() as f32;
            // todo: Consider cacheing the calculations here instead of running this each time.
            for pt in data {
                points.push(
                    to_screen
                        * pos2(
                            pt.0 as f32 / num_pts * width,
                            // Offset for 0 baseline.
                            (pt.1 + MAX_VAL / 2.) / MAX_VAL * CHART_HEIGHT,
                        ),
                );
            }

            let line = Shape::Path(PathShape::line(points, stroke));

            let mut x_axis = Vec::new();
            const NUM_X_TICKS: usize = 12;
            let data_range = data[data.len() - 1].0 - data[0].0;
            let dr_nt = data_range / NUM_X_TICKS;

            let x_axis_posit = CHART_HEIGHT - 4.;

            for i in 0..NUM_X_TICKS {
                let tick_v = data[0].0 + dr_nt * i;

                x_axis.push(ui.ctx().fonts(|fonts| {
                    Shape::text(
                        fonts,
                        to_screen * pos2(i as f32 / NUM_X_TICKS as f32 * width, x_axis_posit),
                        Align2::CENTER_CENTER,
                        tick_v.to_string(),
                        FontId::new(14., FontFamily::Proportional),
                        TICK_COLOR,
                    )
                }));
            }

            ui.painter().extend([line]);
            ui.painter().extend(x_axis);
        });
}

fn draw_proteins(state: &mut State, ui: &mut Ui) {
    for protein in &mut state.volatile.proteins {
        ui.heading(RichText::new(&protein.feature.label).color(Color32::LIGHT_BLUE));

        ui.horizontal(|ui| {
            ui.label(format!(
                "Reading frame: {}, Range: {}",
                protein.reading_frame_match.frame, protein.reading_frame_match.range
            ));
            ui.add_space(COL_SPACING);

            if protein.weight != protein.weight_with_prepost {
                // Only show this segment if pre and post-coding sequences exist
                ui.label(format!(
                    "(Coding region only): AA len: {}  Weight: {:.1}kDa",
                    protein.aa_seq.len(),
                    protein.weight,
                ));
                ui.add_space(COL_SPACING);
            }

            ui.label(format!(
                "AA len: {}  Weight: {:.1}kDa",
                protein.aa_seq.len()
                    + protein.aa_seq_precoding.len()
                    + protein.aa_seq_postcoding.len(),
                protein.weight_with_prepost,
            ));
        });
        ui.add_space(ROW_SPACING / 2.);

        let aa_text = make_aa_text(&protein.aa_seq, state.ui.aa_ident_disp);
        let aa_text_precoding = make_aa_text(&protein.aa_seq_precoding, state.ui.aa_ident_disp);
        let aa_text_postcoding = make_aa_text(&protein.aa_seq_postcoding, state.ui.aa_ident_disp);

        if !aa_text_precoding.is_empty() {
            ui.label(
                RichText::new(aa_text_precoding)
                    .color(COLOR_PRE_POST_CODING_SEQ)
                    .font(FontId::new(FONT_SIZE_SEQ, FontFamily::Monospace)),
            );
        }

        ui.label(
            RichText::new(aa_text)
                .color(COLOR_PROT_SEQ)
                .font(FontId::new(FONT_SIZE_SEQ, FontFamily::Monospace)),
        );

        if !aa_text_postcoding.is_empty() {
            ui.label(
                RichText::new(aa_text_postcoding)
                    .color(COLOR_PRE_POST_CODING_SEQ)
                    .font(FontId::new(FONT_SIZE_SEQ, FontFamily::Monospace)),
            );
        }

        if protein.show_hydropath {
            ui.horizontal(|ui| {
                ui.heading("Hydrophobicity");
                ui.add_space(COL_SPACING);

                if ui.button("Hide hydrophobicity data").clicked() {
                    protein.show_hydropath = false;
                }
            });

            hydrophobicity_chart(&protein.hydropath_data, ui);
        } else {
            if ui.button("Show hydrophobicity data").clicked() {
                protein.show_hydropath = true;
            }
        }

        ui.add_space(ROW_SPACING);
    }
}

pub fn protein_page(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        ui.heading("Proteins, from coding regions");

        ui.add_space(COL_SPACING);
        ui.label("One letter ident:");
        let mut one_letter = state.ui.aa_ident_disp == AaIdent::OneLetter;
        if ui.checkbox(&mut one_letter, "").changed() {
            state.ui.aa_ident_disp = if one_letter {
                AaIdent::OneLetter
            } else {
                AaIdent::ThreeLetters
            };
        }
    });
    ui.add_space(ROW_SPACING);

    if state
        .generic
        .features
        .iter()
        .any(|f| f.feature_type == FeatureType::CodingRegion)
    {
        draw_proteins(state, ui);
    } else {
        ui.label("Create one or more Coding Region feature to display proteins here.");
    }
}
