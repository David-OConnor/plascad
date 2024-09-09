use std::{path::PathBuf, str::FromStr};

use eframe::{
    egui::{pos2, vec2, Frame, Pos2, Rect, RichText, Sense, Shape, Stroke, TextEdit, Ui},
    emath::RectTransform,
    epaint::Color32,
};

use crate::{
    cloning::{make_product_tab, setup_insert_seqs, CloningInsertData},
    file_io::save::load_import,
    gui,
    gui::{
        circle_zoomed::{draw_linear_map, OFFSET},
        navigation::{get_tabs, DEFAULT_TAB_NAME},
        BACKGROUND_COLOR, COL_SPACING, LINEAR_MAP_HEIGHT, ROW_SPACING,
    },
    sequence::{seq_from_str, seq_to_str},
    util::map_linear,
    State,
};

/// Draw a selector for the insert, based on loading from a file.
fn insert_selector(data: &mut CloningInsertData, ui: &mut Ui) {
    for (i, feature) in data.features_loaded.iter().enumerate() {
        let mut border_width = 0.;
        if let Some(j) = data.feature_selected {
            if i == j {
                border_width = 1.;
            }
        }

        Frame::none()
            .stroke(Stroke::new(border_width, Color32::LIGHT_RED))
            .inner_margin(border_width)
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    if ui.button("Select").clicked {
                        data.feature_selected = Some(i);

                        if let Some(seq_this_ft) = feature.range.index_seq(&data.seq_loaded) {
                            seq_this_ft.clone_into(&mut data.seq_insert);
                        }
                    }

                    if !feature.label.is_empty() {
                        ui.label(&feature.label);
                        ui.add_space(COL_SPACING);
                    }

                    let (r, g, b) = feature.feature_type.color();
                    ui.label(
                        RichText::new(feature.feature_type.to_string())
                            .color(Color32::from_rgb(r, g, b)),
                    );
                    ui.add_space(COL_SPACING);

                    ui.label(feature.location_descrip(data.seq_loaded.len()));
                    ui.add_space(COL_SPACING);

                    // +1 because it's inclusive.
                    ui.label(feature.location_descrip(data.seq_loaded.len()));
                });
            });
    }
}

fn insert_file_section(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        ui.label("Choose insert from:");

        // Add buttons for each opened tab
        for (name, i) in get_tabs(
            &state.path_loaded,
            &state.generic[state.active].metadata,
            true,
        ) {
            if ui
                .button(name)
                .on_hover_text("Select an insert from this sequence")
                .clicked()
            {
                let gen = &state.generic[i];
                // This setup, including the break and variables, prevents borrow errors.
                let g = gen.features.clone();
                let s = gen.seq.clone();
                setup_insert_seqs(state, g, s);
                break;
            }
        }

        ui.add_space(COL_SPACING);
        if ui
            .button("Pick insert from file")
            .on_hover_text(
                "Choose a GenBank, PlasCAD, SnapGene, or FASTA file to \
        select an insert from. FASTA files require manual index selection.",
            )
            .clicked()
        {
            state.ui.file_dialogs.cloning_load.select_file();
        }

        ui.add_space(COL_SPACING);

        state.ui.file_dialogs.cloning_load.update(ui.ctx());

        if let Some(path) = state.ui.file_dialogs.cloning_load.take_selected() {
            if let Some(state_loaded) = load_import(&path) {
                // todo: Is there a way to do this without cloning?
                setup_insert_seqs(
                    state,
                    state_loaded.generic.features.clone(),
                    state_loaded.generic.seq.clone(),
                );
            }
        }
    });
}

/// Draw a mini sequence display in its own canvas. Zoomed in on the cloning insert location, and with
/// a line drawn on it.
pub fn seq_lin_disp_cloning(state: &State, ui: &mut Ui) {
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

            const NT_WIDTH: usize = 100;

            // todo: More DRY
            let index_left = (state.cloning_insert_loc as isize - (NT_WIDTH / 2) as isize)
                .rem_euclid(state.get_seq().len() as isize) as usize; // Rust awk % on negative values.
            let index_right = (state.cloning_insert_loc + NT_WIDTH / 2) % state.get_seq().len();

            let mut shapes =
                draw_linear_map(&state, &to_screen, index_left, index_right, false, ui);

            // todo: More DRY.
            let pixel_left = OFFSET.x;
            let pixel_right = ui.available_width() - 2. * OFFSET.x;

            // todo: This segment is DRY from draw_linear_map. Standalone fn etcx. {}    let index_to_x = |mut i: usize| {
            let index_to_x = |mut i: usize| {
                // This handles the case when the zoomed-in view is near the top; the left index
                // will be near the end of the sequence, incorrectly calculating the portion-through in
                // the linear map.
                let right = if index_left > index_right {
                    if i < index_right {
                        // Ie, we are to the right of the origin.
                        i += state.get_seq().len()
                    }

                    index_right + state.get_seq().len()
                } else {
                    index_right
                };

                map_linear(
                    i as f32,
                    (index_left as f32, right as f32),
                    (pixel_left, pixel_right),
                )
            };

            // Draw the insertion site.
            let point_top = pos2(index_to_x(state.cloning_insert_loc), 10.);
            let point_bottom = pos2(index_to_x(state.cloning_insert_loc), 30.);

            shapes.push(Shape::line_segment(
                [to_screen * point_bottom, to_screen * point_top],
                Stroke::new(3., Color32::YELLOW),
            ));
            ui.painter().extend(shapes);
        });
}

pub fn seq_editor_slic(state: &mut State, ui: &mut Ui) {
    ui.heading("SLIC and FastCloning");

    // todo: Once you add this capability.
    ui.label("Clone a sequence into this one. Below, either paste the insert sequence, or select a \
    file (GenBank, SnapGene, FASTA, or PlasCAD) containing the insert sequence. Set the insert location: \
    This is the position in the vector (The currently open sequence) the insert will be placed after. Then click \"Clone\". This will \
    create and open a new sequence, and create optimized primers for both the insert and vector.");

    ui.add_space(ROW_SPACING);

    // ui.label(
    //     "A file dialog will open, prompting you to save a new file for the combined product. Your \
    // current (vector) file will be saved, and the new cloning product file will be opened.",
    // );

    // todo: Zoomed in on insert loc, and draw insert loc.
    seq_lin_disp_cloning(state, ui);
    ui.add_space(ROW_SPACING / 2.);

    ui.horizontal(|ui| {
        ui.label("Insert location: ");
        let mut entry = state.cloning_insert_loc.to_string();
        if ui
            .add(TextEdit::singleline(&mut entry).desired_width(40.))
            .changed()
        {
            state.cloning_insert_loc = entry.parse().unwrap_or(0);
        }

        ui.add_space(COL_SPACING);

        if state.ui.cloning_insert.seq_insert.len() > 6 {
            if ui
                .button(RichText::new("Clone").color(Color32::GOLD))
                .clicked()
            {
                // Save this vector; this file or quicksave instance will be turned into the cloning
                // product.

                make_product_tab(state);
            }
        }
    });

    ui.add_space(ROW_SPACING);
    insert_file_section(state, ui);

    ui.add_space(ROW_SPACING);
    insert_selector(&mut state.ui.cloning_insert, ui);

    ui.add_space(ROW_SPACING);
    ui.horizontal(|ui| {
        ui.heading("Insert:");
        ui.label(&format!(
            "len: {}",
            state.ui.cloning_insert.seq_insert.len()
        ));
    });

    let response = ui.add(
        TextEdit::multiline(&mut state.ui.cloning_insert.seq_input)
            .desired_width(ui.available_width()),
    );
    if response.changed() {
        // Forces only valid NTs to be included in the string.
        state.ui.cloning_insert.seq_insert = seq_from_str(&state.ui.cloning_insert.seq_input);
        state.ui.cloning_insert.seq_input = seq_to_str(&state.ui.cloning_insert.seq_insert);
    }

    ui.add_space(ROW_SPACING);
}
