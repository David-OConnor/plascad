use std::{
    env,
    path::{Path, PathBuf},
    str::FromStr,
};

use eframe::{
    egui::{Frame, RichText, Stroke, TextEdit, Ui},
    epaint::Color32,
};

use crate::{
    file_io::{
        genbank::import_genbank,
        save::{import_fasta, save, StateToSave, DEFAULT_PREFS_FILE, DEFAULT_SAVE_FILE},
        snapgene::import_snapgene,
    },
    gui::{
        navigation::PageSeq, save::save_current_file, set_window_title, COL_SPACING, ROW_SPACING,
    },
    primer::make_cloning_primers,
    sequence::{seq_from_str, seq_to_str, Feature, FeatureDirection, FeatureType},
    util::RangeIncl,
    CloningInsertData, Selection, State,
};
use crate::file_io::save::save_new_product;

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

pub fn seq_editor_slic(state: &mut State, ui: &mut Ui) {
    ui.heading("SLIC and FastCloning");

    // todo: Once you add this capability.
    ui.label("Clone a sequence into this one. Below, either paste the insert sequence, or select a \
    file (GenBank, SnapGene, FASTA, or PlasCAD) containing the insert sequence. Ser the insert location: \
    This is the position in the vector (The currently open sequence) the insert will be placed after. Then click \"Clone\". This will \
    update the sequence with the insert, and create optimized primers for both the insert and vector.");

    ui.add_space(ROW_SPACING);

    ui.label(
        "A file dialog will open, prompting you to save a new file for the combined product. Your \
    current (vector) file will be saved, and the new cloning product file will be opened.",
    );

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
            // state.ui.file_dialogs.selected = Some(path.to_owned());

            if let Some(extension) = path.extension().and_then(|ext| ext.to_str()) {
                match extension.to_lowercase().as_ref() {
                    "pcad" => {
                        let data =
                            State::load(&path, &PathBuf::from_str(DEFAULT_PREFS_FILE).unwrap());
                        state.ui.cloning_insert.features_loaded = data.generic.features;
                        state.ui.cloning_insert.seq_loaded = data.generic.seq;
                    }
                    "fasta" => {
                        if let Ok((seq, _id, _description)) = import_fasta(&path) {
                            state.ui.cloning_insert.seq_loaded = seq;
                        }
                    }
                    "dna" => {
                        if let Ok(data) = import_snapgene(&path) {
                            state.ui.cloning_insert.features_loaded = data.features;
                            state.ui.cloning_insert.seq_loaded = data.seq;
                        }
                    }
                    "gb" | "gbk" => {
                        if let Ok(data) = import_genbank(&path) {
                            state.ui.cloning_insert.features_loaded = data.features;
                            state.ui.cloning_insert.seq_loaded = data.seq;
                        }
                    }
                    _ => {
                        eprintln!(
                            "The file to import must be in FASTA, GenBank, or SnapGene format."
                        )
                    }
                }

                // Choose the initial insert as the CDS or gene with the largest len.
                let mut best = None;
                let mut best_len = 0;
                for (i, feature) in state.ui.cloning_insert.features_loaded.iter().enumerate() {
                    let len = feature.len(state.generic.seq.len());
                    if (feature.feature_type == FeatureType::CodingRegion
                        || feature.feature_type == FeatureType::Gene)
                        && len > best_len
                    {
                        best_len = len;
                        best = Some(i);
                    }
                }

                if let Some(feat_i) = best {
                    let feature = &state.ui.cloning_insert.features_loaded[feat_i];

                    if let Some(seq_this_ft) =
                        feature.range.index_seq(&state.ui.cloning_insert.seq_loaded)
                    {
                        state.ui.cloning_insert.feature_selected = best;
                        state.ui.cloning_insert.seq_insert = seq_this_ft.to_owned();
                        state.ui.cloning_insert.seq_input = seq_to_str(&seq_this_ft);
                    }
                }
            }
        }

        if state.ui.cloning_insert.seq_insert.len() > 6 {
            if ui.button(RichText::new("Clone").color(Color32::GOLD)).clicked() {
                // Save this vector; this file or quicksave instance will be turned into the cloning
                // product.
                save_current_file(state);

                // Make sure to create cloning primers before performing the insert, or the result will be wrong.
                make_cloning_primers(state);

                // todo: Unecessary clone
                let insert = state.ui.cloning_insert.seq_insert.clone();
                state.insert_nucleotides(&insert, state.cloning_insert_loc);

                let label = match state.ui.cloning_insert.feature_selected {
                    Some(i) => state.ui.cloning_insert.features_loaded[i].label.clone(),
                    None => "Cloning insert".to_owned(),
                };

                // todo: Eventually, we'll likely be pulling in sequences already associated with a feature;
                // todo: Use the already existing data instead.
                state.generic.features.push(Feature {
                    range: RangeIncl::new(
                        state.cloning_insert_loc + 1,
                        state.cloning_insert_loc + state.ui.cloning_insert.seq_insert.len(),
                    ),
                    label,
                    feature_type: FeatureType::CodingRegion,
                    direction: FeatureDirection::Forward,
                    ..Default::default()
                });

                state.ui.page_seq = PageSeq::View;
                state.ui.selected_item = Selection::Feature(state.generic.features.len() - 1);

                save_new_product("SLIC cloning product", state, ui);
            }
        }
    });

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
