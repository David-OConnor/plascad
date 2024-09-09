//! GUI code related to ligation operations.

use eframe::{
    egui::{
        pos2, vec2, Align2, Color32, FontFamily, FontId, Frame, Pos2, Rect, RichText, ScrollArea,
        Sense, Shape, Stroke, Ui,
    },
    emath::RectTransform,
};

use crate::{
    gui::{
        circle::{FEATURE_OUTLINE_COLOR, FEATURE_STROKE_WIDTH},
        navigation::get_tabs,
        seq_lin_disp, BACKGROUND_COLOR, COL_SPACING, ROW_SPACING,
    },
    ligation::{digest, ligate, LigationFragment},
    sequence::seq_to_str,
    util::map_linear,
    State,
};

// This X offset must have room for the RE Nts displayed on the left.
const OFFSET_X: f32 = 80.;
const OFFSET_Y: f32 = 30.;

const PRODUCT_ROW_SPACING: f32 = 60.;
const FRAG_DISP_HEIGHT: f32 = 30.;
const FRAG_DISP_HEIGHT_DIV2: f32 = FRAG_DISP_HEIGHT / 2.;

/// Draw a graphical depiction of digestion products.
fn draw_graphics(products: &[LigationFragment], seq_len: usize, ui: &mut Ui) {
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

            // Scale pixel length based on the full sequence length.
            let left_edge_px = OFFSET_X;
            let right_edge_px = ui.available_width() - OFFSET_X - 20.;

            for frag in products {
                let start_x = OFFSET_X;
                let end_x = start_x + frag.seq.len() as f32 * 2.; // todo: Rough.

                // todo: This approach won't work when adding in fragments from other sequences.
                let end_x = map_linear(
                    frag.seq.len() as f32,
                    (0., seq_len as f32),
                    (left_edge_px, right_edge_px),
                );

                // println!("VEC: {:?}",                    vec![
                //     pos2(start_x, row_px - FRAG_DISP_HEIGHT_DIV2),
                //     pos2(end_x, row_px - FRAG_DISP_HEIGHT_DIV2),
                //     pos2(end_x, row_px + FRAG_DISP_HEIGHT_DIV2),
                //     pos2(start_x, row_px + FRAG_DISP_HEIGHT_DIV2),
                // ]);

                // return; // todo temp

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

                let box_center = pos2((end_x + start_x) / 2., row_px);

                // Draw the RE ends.
                let label_pt_left_top = pos2(start_x - 10., row_px - 10.);
                let label_pt_right_top = pos2(end_x + 10., row_px - 10.);

                let label_pt_left_bottom = pos2(start_x - 10., row_px + 10.);
                let label_pt_right_bottom = pos2(end_x + 10., row_px + 10.);

                let (re_text_left_top, re_text_left_bottom, re_name_left) = match &frag.re_left {
                    Some(re) => (
                        seq_to_str(&re.overhang_top_left()),
                        seq_to_str(&re.overhang_bottom_left()),
                        re.name.clone(),
                    ),
                    None => (String::new(), String::new(), String::new()),
                };

                let (re_text_right_top, re_text_right_bottom, re_name_right) = match &frag.re_right
                {
                    Some(re) => (
                        seq_to_str(&re.overhang_top_right()),
                        seq_to_str(&re.overhang_bottom_right()),
                        re.name.clone(),
                    ),
                    None => (String::new(), String::new(), String::new()),
                };

                // todo: Allow viewing and copying fragment seqs, eg from selecting them.

                // Draw info like fragment length inside the boxes.
                shapes.push(ui.ctx().fonts(|fonts| {
                    Shape::text(
                        fonts,
                        to_screen * box_center,
                        Align2::CENTER_CENTER,
                        &format!("{} bp", frag.seq.len()),
                        FontId::new(16., FontFamily::Proportional),
                        Color32::BLACK,
                    )
                }));

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

                // Text description of which enzymes were used.
                let enzyme_text = format!("{}  |  {}", re_name_left, re_name_right);
                let enzyme_text_pos = to_screen * pos2(end_x + 70., row_px);

                shapes.push(ui.ctx().fonts(|fonts| {
                    Shape::text(
                        fonts,
                        enzyme_text_pos,
                        Align2::LEFT_CENTER,
                        &enzyme_text,
                        FontId::new(16., FontFamily::Proportional),
                        Color32::LIGHT_RED,
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
        ui.heading("Digestion and ligation");
        ui.add_space(COL_SPACING * 2.);

        ui.label("Unique cutters only:");
        ui.checkbox(&mut state.ui.re.unique_cutters_only, "");
        ui.add_space(COL_SPACING * 2.);

        ui.label("Sticky ends only:");
        ui.checkbox(&mut state.ui.re.sticky_ends_only, "");
    });

    seq_lin_disp(state, ui, true);

    ui.add_space(ROW_SPACING);

    // let mut res_matched = HashMap::new(); // Re, match count
    let mut res_matched = Vec::new();
    for re_match in &state.volatile[state.active].restriction_enzyme_matches {
        if re_match.lib_index + 1 >= state.restriction_enzyme_lib.len() {
            continue;
        }

        let re = &state.restriction_enzyme_lib[re_match.lib_index];

        if !res_matched.contains(&re) {
            if (state.ui.re.unique_cutters_only && re_match.match_count > 1)
                || (state.ui.re.sticky_ends_only && re.makes_blunt_ends())
            {
                continue;
            }
            res_matched.push(re);
        }
    }

    ui.horizontal(|ui| {
        ui.heading("Opened files to digest");
        ui.add_space(COL_SPACING);

        for (name, i) in get_tabs(
            &state.path_loaded,
            &state.generic[state.active].metadata,
            true,
        ) {
            // todo: DRY with page selectors and Cloning.
            let color = if state.ui.re.tabs_selected.contains(&i) {
                Color32::GREEN
            } else {
                Color32::WHITE
            };
            let button = ui.button(
                RichText::new(name).color(color), // .background_color(TAB_BUTTON_COLOR),
            );

            if button.clicked() {
                if state.ui.re.tabs_selected.contains(&i) {
                    // todo: DRY with res_selected below.
                    for (j, tab_i) in state.ui.re.tabs_selected.iter().enumerate() {
                        if i == *tab_i {
                            state.ui.re.tabs_selected.remove(j);
                            break;
                        }
                    }
                } else {
                    state.ui.re.tabs_selected.push(i);
                }
            }

            ui.add_space(COL_SPACING / 2.);
        }
    });
    ui.add_space(ROW_SPACING);

    // todo: Highlight common (and later, compatible) RE matches among fragments.

    ui.horizontal(|ui| {
        ui.heading("Restriction enzymes matched");
        ui.add_space(COL_SPACING);

        ui.label("Click to select for digestion.");
    });

    let res_per_row = 8; // todo: Based on screen width etc.
    let num_re_rows = res_matched.len() / res_per_row + 1;

    for row_i in 0..num_re_rows {
        // todo: YOu must make this wrap, or take a different approach.
        ui.horizontal(|ui| {
            for col_i in 0..res_per_row {
                let index = row_i * res_per_row + col_i;
                if index + 1 > res_matched.len() {
                    break;
                }
                let re = res_matched[index];

                let selected = state.ui.re.res_selected.contains(&re.name);

                let color = if selected {
                    Color32::LIGHT_GREEN
                } else {
                    Color32::WHITE
                };

                if ui.button(RichText::new(&re.name).color(color)).clicked {
                    if selected {
                        for (i, name) in state.ui.re.res_selected.iter().enumerate() {
                            if name == &re.name {
                                state.ui.re.res_selected.remove(i);
                                break;
                            }
                        }
                    } else {
                        state.ui.re.res_selected.push(re.name.clone());
                    }
                }

                ui.add_space(COL_SPACING / 2.);
                ui.label(re.cut_depiction());

                //
                // if !state.ui.re.unique_cutters_only { // No point in displaying count for unique cutters; always 1.
                //      // todo: Show count.
                //      // ui.label(format!(": {}"));
                // }
            }
            ui.add_space(COL_SPACING);
        });
    }

    ui.add_space(ROW_SPACING);

    ui.horizontal(|ui| {
        if !state.ui.re.res_selected.is_empty() {
            if ui
                .button(RichText::new("Digest").color(Color32::GOLD))
                .clicked()
            {
                state.volatile[state.active].re_digestion_products = digest(
                    &state.ui.re.res_selected,
                    &state.volatile[state.active].restriction_enzyme_matches,
                    &state.restriction_enzyme_lib,
                    state.get_seq(),
                    state.generic[state.active].topology,
                );
            }
        }

        if !state.volatile[state.active]
            .re_digestion_products
            .is_empty()
        {
            ui.add_space(COL_SPACING);
            if ui
                .button(RichText::new("Ligate").color(Color32::GOLD))
                .clicked()
            {
                // state.volatile[state.active].re_ligation_products = ligate(
                //     &state.volatile[state.active].re_digestion_products,
                // );
            }
        }
    });

    // todo: YOu need to draw complementary strands and overhangs.
    // todo: Likely just rectangles for the strands, with the RE sites as the only letters. For both ends.

    // for frag in &state.volatile[state.active].re_digestion_products {
    //     ui.label(format!("{} - {}", frag.seq.len(), seq_to_str(&frag.seq)));
    // }

    ScrollArea::vertical().id_source(100).show(ui, |ui| {
        draw_graphics(
            &state.volatile[state.active].re_digestion_products,
            state.get_seq().len(),
            ui,
        );
    });
}
