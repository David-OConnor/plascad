use std::{path::PathBuf, str::FromStr};

use eframe::{
    egui::{Frame, RichText, Stroke, TextEdit, Ui},
    epaint::Color32,
};

use crate::{
    cloning::{make_product_tab, setup_insert_seqs, CloningInsertData},
    file_io::save::load_import,
    gui::{navigation::DEFAULT_TAB_NAME, COL_SPACING, ROW_SPACING},
    sequence::{seq_from_str, seq_to_str},
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
                            data.seq_insert = seq_this_ft.to_owned();
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
        for gen in &mut state.generic {
            let name = &gen.metadata.plasmid_name;
            let btn_name = if name.len() > 20 {
                format!("{}...", &name[..20].to_string())
            } else if name.is_empty() {
                DEFAULT_TAB_NAME.to_owned()
            } else {
                name.to_owned()
            };

            if ui
                // todo: Consider something like the actual tab name system, that uses file name,
                // todo or a default A/R
                .button(btn_name)
                .on_hover_text("Select an insert from this sequence")
                .clicked()
            {
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

pub fn seq_editor_slic(state: &mut State, ui: &mut Ui) {
    ui.heading("SLIC and FastCloning");

    // todo: Once you add this capability.
    ui.label("Clone a sequence into this one. Below, either paste the insert sequence, or select a \
    file (GenBank, SnapGene, FASTA, or PlasCAD) containing the insert sequence. Ser the insert location: \
    This is the position in the vector (The currently open sequence) the insert will be placed after. Then click \"Clone\". This will \
    update the sequence with the insert, and create optimized primers for both the insert and vector.");

    ui.add_space(ROW_SPACING);

    // ui.label(
    //     "A file dialog will open, prompting you to save a new file for the combined product. Your \
    // current (vector) file will be saved, and the new cloning product file will be opened.",
    // );

    ui.add_space(ROW_SPACING);

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
