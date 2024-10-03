//! GUI code related to ligation operations.

use std::{collections::HashMap, path::PathBuf};

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
        lin_maps::seq_lin_disp,
        navigation::get_tabs,
        select_color_text,
        theme::COLOR_ACTION,
        BACKGROUND_COLOR, COL_SPACING, ROW_SPACING,
    },
    ligation,
    ligation::{digest, LigationFragment},
    sequence::seq_to_str,
    util::{map_linear, name_from_path},
    State,
};

// This X offset must have room for the RE Nts displayed on the left.
const OFFSET_X: f32 = 60.;
const OFFSET_Y: f32 = 30.;

const PRODUCT_ROW_SPACING: f32 = 60.;
const FRAG_DISP_HEIGHT: f32 = 30.;
const FRAG_DISP_HEIGHT_DIV2: f32 = FRAG_DISP_HEIGHT / 2.;

// We cycle through this for visual separation
const SEQ_COLORS: [Color32; 4] = [
    Color32::LIGHT_RED,
    Color32::LIGHT_GREEN,
    Color32::LIGHT_YELLOW,
    Color32::LIGHT_GRAY,
];

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

            // let rect_size = response.rect.size();

            let mut shapes = Vec::new();

            let stroke = Stroke::new(FEATURE_STROKE_WIDTH, FEATURE_OUTLINE_COLOR);

            let mut row_px = OFFSET_Y;

            // Scale pixel length based on the full sequence length.
            let left_edge_px = OFFSET_X;
            // Leave enough room on both sides for descriptive text annotations.
            let right_edge_px = ui.available_width() - OFFSET_X - 80.;

            let mut seq_colors: HashMap<String, Color32> = HashMap::new();
            let mut color_i = 0;

            for frag in products {
                let frag_color = if seq_colors.contains_key(&frag.source_name) {
                    seq_colors[&frag.source_name]
                } else {
                    color_i = (color_i + 1) % SEQ_COLORS.len();
                    seq_colors.insert(frag.source_name.to_owned(), SEQ_COLORS[color_i]);
                    SEQ_COLORS[color_i]
                };

                let start_x = OFFSET_X;
                // let end_x = start_x + frag.seq.len() as f32 * 2.; // todo: Rough.

                // todo: This approach won't work when adding in fragments from other sequences.
                let end_x = map_linear(
                    frag.seq.len() as f32,
                    (0., seq_len as f32),
                    (left_edge_px, right_edge_px),
                );

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
                        FontId::new(14., FontFamily::Proportional),
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
                let enzyme_text_pos = to_screen * pos2(end_x + 50., row_px);
                let source_name_text_pos = to_screen * pos2(end_x + 150., row_px);

                shapes.push(ui.ctx().fonts(|fonts| {
                    Shape::text(
                        fonts,
                        enzyme_text_pos,
                        Align2::LEFT_CENTER,
                        &enzyme_text,
                        FontId::new(14., FontFamily::Proportional),
                        Color32::LIGHT_RED,
                    )
                }));

                // Show the source name this fragment came from.

                shapes.push(ui.ctx().fonts(|fonts| {
                    Shape::text(
                        fonts,
                        source_name_text_pos,
                        Align2::LEFT_CENTER,
                        &frag.source_name,
                        FontId::new(14., FontFamily::Proportional),
                        frag_color,
                    )
                }));

                row_px += PRODUCT_ROW_SPACING;
            }

            ui.painter().extend(shapes);
        });
}

fn tab_selection(
    tabs: &mut Vec<usize>,
    path_loaded: &[Option<PathBuf>],
    plasmid_names: &[&str],
    ui: &mut Ui,
) -> bool {
    let mut clear_res = false;

    ui.horizontal(|ui| {
        ui.heading("Opened files to digest");
        ui.add_space(COL_SPACING);

        for (name, i) in get_tabs(path_loaded, plasmid_names, true) {
            // todo: DRY with page selectors and Cloning.
            let selected = tabs.contains(&i);

            let button = ui.button(
                select_color_text(&name, selected), // .background_color(TAB_BUTTON_COLOR),
            );

            if button.clicked() {
                if tabs.contains(&i) {
                    // todo: DRY with res_selected below.
                    for (j, tab_i) in tabs.iter().enumerate() {
                        if i == *tab_i {
                            tabs.remove(j);
                            // This is crude; we only need to re-select REs that were enabled by this tab.
                            clear_res = true;
                            break;
                        }
                    }
                } else {
                    tabs.push(i);
                }
            }

            ui.add_space(COL_SPACING / 2.);
        }
    });
    clear_res
}

// // todo
// /// Remove selected res if they are filtered out.
// fn filter_res_selected(data: &mut ReUi, lib: &[RestrictionEnzyme]) {
//     for name in &mut data.res_selected {
//         let mut re = None;
//         for re_ in lib {
//             // if &re_.name == name {
//             if &re_.name == name {
//                 re = Some(re_.clone());
//                 break;
//             }
//             if re.is_none() {
//                 eprintln!("Error: Unable to find matching RE");
//                 break;
//             }
//         }
//
//         let re = re.unwrap();
//     }
// }

pub fn ligation_page(state: &mut State, ui: &mut Ui) {
    // todo: Scrolling is not working
    ScrollArea::vertical().id_source(100).show(ui, |ui| {
        // todo: Cache calcs in this fn A/R
        // Adjust the nucleotide width of digest bars based on the longest sequence selected.
        let mut longest_seq = 0;
        for active in &state.ui.re.tabs_selected {
            if state.generic[*active].seq.len() > longest_seq {
                longest_seq = state.generic[*active].seq.len();
            }
        }

        // todo: Don't run this every frame! Store in state, and update it only when something notable changes.
        let res_matched = ligation::filter_res(&state.ui.re, &state.volatile, &state.restriction_enzyme_lib);

        ui.horizontal(|ui| {
            ui.heading("Digestion and ligation");
            ui.add_space(COL_SPACING * 2.);

            let mut changed_filters = false;

            ui.label("Unique cutters only:").on_hover_text("Only display restriction enzymes that cut a sequence exactly once. (Affects display on other pages as well).");
            if ui.checkbox(&mut state.ui.re.unique_cutters_only, "").changed() {
                changed_filters = true;
            };

            ui.add_space(COL_SPACING);

            ui.label("Sticky ends only:").on_hover_text("Only display restriction enzymes that produce overhangs (sticky ends). Don't show ones that produce blunt ends. (Affects display on other pages as well).");
            if ui.checkbox(&mut state.ui.re.sticky_ends_only, "").changed() {
                changed_filters = true;
            };

            ui.add_space(COL_SPACING);

            ui.label("2+ sequences only:");
            if ui.checkbox(&mut state.ui.re.multiple_seqs, "")
                .on_hover_text("Only display restriction enzymes that cut two or more selected sequences, if applicable.")
                .changed() {
                changed_filters = true;
            };

            ui.add_space(COL_SPACING);

            // If we've changed filters, update REs selected IRT these filtersx.
            if changed_filters {
                // filter_res_selected(&mut state.ui.re, &state.restriction_enzyme_lib);
            }

        });
        ui.add_space(ROW_SPACING);

        let plasmid_names: &Vec<_> = &state.generic.iter().map(|v| v.metadata.plasmid_name.as_str()).collect();

        let clear_res = tab_selection(&mut state.ui.re.tabs_selected, &state.path_loaded, plasmid_names, ui);
        if clear_res {
            state.ui.re.res_selected = Vec::new();
        }

        ui.add_space(ROW_SPACING / 2.);

        // // todo: This is potentially un-ideal computationally. Note that the name-based code is stored permanently.
        // // todo: NOte that this is just for the highlights. todo: Just store teh whole Re Match instead of names as state.ui.res_selected.
        // let mut res_sel = Vec::new();
        // for name in &state.ui.re.res_selected {
        //     for re in &state.restriction_enzyme_lib {
        //         if name == &re.name {
        //             res_sel.push(re);
        //         }
        //     }
        // }

        for active in &state.ui.re.tabs_selected {
            seq_lin_disp(&state.generic[*active], true, state.ui.selected_item, &state.ui.re.res_selected, None, &state.ui,&state.volatile[state.active].restriction_enzyme_matches,
                         &state.restriction_enzyme_lib, ui);
            ui.add_space(ROW_SPACING/2.);
        }

        // todo: Highlight common (and later, compatible) RE matches among fragments.

        ui.horizontal(|ui| {
            ui.heading("Restriction enzymes matched");
            ui.add_space(COL_SPACING);

            ui.label("Click to select for digestion.");
        });

        let res_per_row = 8; // todo: Based on screen width etc.
        let num_re_rows = res_matched.len() / res_per_row + 1;

        for row_i in 0..num_re_rows {
            ui.horizontal(|ui| {
                for col_i in 0..res_per_row {
                    let index = row_i * res_per_row + col_i;
                    if index + 1 > res_matched.len() {
                        break;
                    }
                    let re = res_matched[index];

                    let selected = state.ui.re.res_selected.contains(&re);

                    if ui.button(select_color_text(&re.name, selected)).clicked {
                        if selected {
                            for (i, re_sel) in state.ui.re.res_selected.iter().enumerate() {
                                if re_sel == re {
                                    state.ui.re.res_selected.remove(i);
                                    break;
                                }
                            }
                        } else {
                            state.ui.re.res_selected.push(re.clone());
                        }
                    }
                    ui.label(re.cut_depiction());
                    ui.add_space(COL_SPACING / 2.);
                }
                ui.add_space(COL_SPACING);
            });
        }

        ui.add_space(ROW_SPACING);

        ui.horizontal(|ui| {
            if !state.ui.re.res_selected.is_empty() {
                if ui
                    .button(RichText::new("Digest").color(COLOR_ACTION))
                    .clicked()
                {
                    let mut products = Vec::new();

                    for active in &state.ui.re.tabs_selected {
                        let source_name = name_from_path(
                            &state.path_loaded[*active],
                            &state.generic[*active].metadata.plasmid_name,
                            true,
                        );
                        products.extend(digest(
                            &source_name,
                            &state.ui.re.res_selected,
                            &state.volatile[*active].restriction_enzyme_matches,
                            &state.restriction_enzyme_lib,
                            &state.generic[*active].seq,
                            state.generic[*active].topology,
                        ));
                    }

                    state.volatile[state.active].re_digestion_products = products;
                }
            }

            if !state.volatile[state.active]
                .re_digestion_products
                .is_empty()
            {
                // todo: Put back when ready
                // ui.add_space(COL_SPACING);
                // if ui
                //     .button(RichText::new("Ligate").color(COLOR_ACTION))
                //     .clicked()
                // {
                //     // state.volatile[state.active].re_ligation_products = ligate(
                //     //     &state.volatile[state.active].re_digestion_products,
                //     // );
                // }
            }
        });

        // Display the digestion products,
        draw_graphics(
            &state.volatile[state.active].re_digestion_products,
            longest_seq,
            ui,
        );
    });
}
