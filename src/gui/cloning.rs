use std::{path::PathBuf, str::FromStr};

use eframe::egui::{RichText, TextEdit, Ui};
use eframe::epaint::Color32;
use crate::{
    file_io::{
        genbank::import_genbank,
        save::{import_fasta, DEFAULT_PREFS_FILE},
        snapgene::import_snapgene,
    },
    gui::{navigation::PageSeq, COL_SPACING, ROW_SPACING},
    primer::make_cloning_primers,
    sequence::{seq_from_str, seq_to_str, Feature, FeatureType},
    Selection, State,
};

pub fn seq_editor_slic(state: &mut State, ui: &mut Ui) {
    ui.heading("SLIC and FastCloning");

    // todo: Once you add this capability.
    ui.label("Clone a sequence into this one. Below, either paste the insert sequence, or select a \
    file (GenBank, SnapGene, FASTA, or PlasCAD) containing the insert sequence. Select the insert location: \
    This is the position in the vector the insert will be placed after. Then click \"Make cloning primers\". This will \
    update the sequence with the insert, and create optimized primers for both the insert and vector.");

    ui.add_space(ROW_SPACING);

    ui.horizontal(|ui| {
        ui.label("Insert location: ");
        let mut entry = state.insert_loc.to_string();
        if ui
            .add(TextEdit::singleline(&mut entry).desired_width(40.))
            .changed()
        {
            state.insert_loc = entry.parse().unwrap_or(0);
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
            state.ui.file_dialogs.selected = Some(path.to_owned());

            if let Some(extension) = path.extension().and_then(|ext| ext.to_str()) {
                match extension.to_lowercase().as_ref() {
                    "pcad" => {
                        let data =
                            State::load(&path, &PathBuf::from_str(DEFAULT_PREFS_FILE).unwrap());
                        state.ui.cloning_insert.features_loaded = data.generic.features;
                        state.ui.cloning_insert.seq_loaded = data.generic.seq;
                    }
                    "fasta" => {
                        if let Ok((seq, id, description)) = import_fasta(&path) {
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
            }
        }

        if ui.button("➕ Make cloning primers").clicked() {
            state.sync_cloning_product();
            state.sync_seq_related(None);
            make_cloning_primers(state);

            // todo: Eventually, we'll likely be pulling in sequences already associated with a feature;
            // todo: Use the already existing data instead.
            state.generic.features.push(Feature {
                index_range: (
                    state.insert_loc + 1,
                    state.insert_loc + 1 + state.ui.cloning_seq_insert_input.len(),
                ), // todo: Check off-by-one.
                label: "Cloning insert".to_owned(),
                feature_type: FeatureType::CodingRegion,
                ..Default::default()
            });

            state.ui.page_seq = PageSeq::View;
            state.ui.selected_item = Selection::Feature(state.generic.features.len() - 1);
        }
    });

    // todo: Select red as you do on teh table.
    for (i, feature) in state.ui.cloning_insert.features_loaded.iter().enumerate() {
        ui.horizontal(|ui| {
            if ui.button("Select").clicked {
                state.ui.cloning_insert.feature_selected = Some(i);
            }

            if !feature.label.is_empty() {
                ui.label(&feature.label);
                ui.add_space(COL_SPACING);
            }

            let (r, g, b) = feature.feature_type.color();
            ui.label(RichText::new(&feature.feature_type.to_string()).color(Color32::from_rgb(r, g, b)));
            ui.add_space(COL_SPACING);

            ui.label(format!("{}..{}", feature.index_range.0, feature.index_range.1));
        });

    }

    ui.horizontal(|ui| {
        ui.heading("Insert:");
        ui.label(&format!("len: {}", state.ui.cloning_seq_insert_input.len()));
    });

    let response = ui.add(
        TextEdit::multiline(&mut state.ui.cloning_seq_insert_input)
            .desired_width(ui.available_width()),
    );
    if response.changed() {
        let seq = seq_from_str(&state.ui.cloning_seq_insert_input);
        state.ui.cloning_seq_insert_input = seq_to_str(&seq);
    }

    ui.add_space(ROW_SPACING);

}
