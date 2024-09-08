//! GUI code related to ligation operations.

use eframe::{
    egui::{pos2, vec2, Color32, Frame, Pos2, Rect, RichText, Sense, Shape, Stroke, Ui},
    emath::RectTransform,
};
use eframe::egui::{Align2, FontFamily, FontId};
use crate::{
    gui::{
        circle::{FEATURE_OUTLINE_COLOR, FEATURE_STROKE_WIDTH},
        BACKGROUND_COLOR, COL_SPACING, ROW_SPACING,
    },
    ligation::digest,
    sequence::seq_to_str,
    State,
};
use crate::ligation::LigationFragment;

// This X offset must have room for the RE Nts displayed on the left.
const OFFSET_X: f32 = 80.;
const OFFSET_Y: f32 = 30.;

const PRODUCT_ROW_SPACING: f32 = 60.;
const FRAG_DISP_HEIGHT: f32 = 30.;
const FRAG_DISP_HEIGHT_DIV2: f32 = FRAG_DISP_HEIGHT / 2.;

/// Draw a graphical depiction of digestion products.
fn draw_graphics(products: &[LigationFragment], ui: &mut Ui) {
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

            let mut shapes = Vec::new();

            let frag_color = Color32::LIGHT_BLUE; // todo A/R

            let stroke = Stroke::new(FEATURE_STROKE_WIDTH, FEATURE_OUTLINE_COLOR);

            let mut row_px = OFFSET_Y;


            for frag in products {
                let start_x = OFFSET_X;
                let end_x = frag.seq.len() as f32 * 2.; // todo: Rough.

                shapes.push(Shape::convex_polygon(
                    vec![
                        to_screen * pos2(start_x, row_px - FRAG_DISP_HEIGHT_DIV2),
                        to_screen * pos2(end_x, row_px - FRAG_DISP_HEIGHT_DIV2),
                        to_screen * pos2(end_x, row_px + FRAG_DISP_HEIGHT_DIV2),
                        to_screen * pos2(start_x, row_px + FRAG_DISP_HEIGHT_DIV2),
                    ],
                    frag_color,
                    stroke,
                ));


                // Draw the RE ends.
                let label_pt_left_top = pos2(start_x - 10., row_px - 10.);
                let label_pt_right_top = pos2(end_x + 10., row_px - 10.);

                let label_pt_left_bottom = pos2(start_x - 10., row_px + 10.);
                let label_pt_right_bottom = pos2(end_x + 10., row_px + 10.);

                let re_text_left_top = match &frag.re_left {
                    Some(re) => seq_to_str(&re.overhang_top_left()),
                    None => String::new()
                };

                let re_text_right_top = match &frag.re_right {
                    Some(re) => seq_to_str(&re.overhang_top_right()),
                    None => String::new()
                };

                let re_text_left_bottom = match &frag.re_left {
                    Some(re) => seq_to_str(&re.overhang_bottom_left()),
                    None => String::new()
                };

                let re_text_right_bottom = match &frag.re_right {
                    Some(re) => seq_to_str(&re.overhang_bottom_right()),
                    None => String::new()
                };

                // Left, top
                shapes.push(ui.ctx().fonts(|fonts| {
                    Shape::text(
                        fonts,
                        to_screen * label_pt_left_top,
                        Align2::RIGHT_CENTER,
                        &re_text_left_top,
                        FontId::new(14., FontFamily::Proportional),
                        Color32::LIGHT_YELLOW,
                    )
                }));

                // Right, top
                shapes.push(ui.ctx().fonts(|fonts| {
                    Shape::text(
                        fonts,
                        to_screen * label_pt_right_top,
                        Align2::LEFT_CENTER,
                        &re_text_right_top,
                        FontId::new(14., FontFamily::Proportional),
                        Color32::LIGHT_YELLOW,
                    )
                }));

                // Left, bottom
                shapes.push(ui.ctx().fonts(|fonts| {
                    Shape::text(
                        fonts,
                        to_screen * label_pt_left_bottom,
                        Align2::RIGHT_CENTER,
                        &re_text_left_bottom,
                        FontId::new(14., FontFamily::Proportional),
                        Color32::LIGHT_YELLOW,
                    )
                }));

                // Right, bottom
                shapes.push(ui.ctx().fonts(|fonts| {
                    Shape::text(
                        fonts,
                        to_screen * label_pt_right_bottom,
                        Align2::LEFT_CENTER,
                        &re_text_right_bottom,
                        FontId::new(14., FontFamily::Proportional),
                        Color32::LIGHT_YELLOW,
                    )
                }));

                row_px += PRODUCT_ROW_SPACING;
            }

            ui.painter().extend(shapes);
        });
}

pub fn ligation_page(state: &mut State, ui: &mut Ui) {
    // todo: Cache this A/R
    ui.horizontal(|ui| {
        ui.heading("Digestion and ligration");
        ui.add_space(COL_SPACING * 2.);

        ui.label("Unique cutters only:");
        ui.checkbox(&mut state.ui.re.unique_cutters_only, "");
    });

    ui.add_space(ROW_SPACING);

    // todo: Only show unique cutters A/R. And/or in the list of sites.
    // let mut res_matched = HashMap::new(); // Re, match count
    let mut res_matched = Vec::new();
    for re_match in &state.volatile[state.active].restriction_enzyme_matches {
        if re_match.lib_index + 1 >= state.restriction_enzyme_lib.len() {
            continue;
        }

        let re = &state.restriction_enzyme_lib[re_match.lib_index];

        if !res_matched.contains(&re) {
            if state.ui.re.unique_cutters_only && re_match.match_count > 1 {
                continue;
            }
            res_matched.push(re);
        }
    }

    ui.horizontal(|ui| {
        ui.heading("Restriction enzymes matched");
        ui.add_space(COL_SPACING);

        ui.label("Click to select for digestion.");
    });

    for re in res_matched {
        let selected = state.ui.re.selected.contains(&re.name);

        let color = if selected {
            Color32::LIGHT_GREEN
        } else {
            Color32::WHITE
        };

        ui.horizontal(|ui| {
            if ui.button(RichText::new(&re.name).color(color)).clicked {
                if selected {
                    for (i, name) in state.ui.re.selected.iter().enumerate() {
                        if name == &re.name {
                            state.ui.re.selected.remove(i);
                            break;
                        }
                    }
                } else {
                    state.ui.re.selected.push(re.name.clone());
                }
            }

            ui.add_space(COL_SPACING);
            ui.label(format!("{}", re.cut_depiction()));
        });

        if !state.ui.re.unique_cutters_only { // No point in displaying count for unique cutters; always 1.
            // todo: Show count.
            // ui.label(format!(": {}"));
        }
    }

    if !state.ui.re.selected.is_empty() {
        ui.add_space(ROW_SPACING);
        if ui
            .button(RichText::new("Digest").color(Color32::GOLD))
            .clicked()
        {
            state.volatile[state.active].re_digestion_products = digest(
                &state.ui.re.selected,
                &state.volatile[state.active].restriction_enzyme_matches,
                &state.restriction_enzyme_lib,
                state.get_seq(),
                state.generic[state.active].topology,
            );
        }
    }

    // todo: YOu need to draw complementary strands and overhangs.
    // todo: Likely just rectangles for the strands, with the RE sites as the only letters. For both ends.

    for frag in &state.volatile[state.active].re_digestion_products {
        ui.label(format!("{} - {}", frag.seq.len(), seq_to_str(&frag.seq)));
    }

    draw_graphics(&state.volatile[state.active].re_digestion_products , ui);
}
