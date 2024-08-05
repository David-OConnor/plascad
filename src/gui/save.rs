//! GUI code for saving and loading

use std::{
    env, io,
    path::{Path, PathBuf},
};

use eframe::egui::Ui;
// use egui_file::FileDialog;
use egui_file_dialog::FileDialog;

use crate::{
    genbank_parse::import_genbank,
    primer::PrimerData,
    save::{export_fasta, import_fasta, load, save, StateToSave, DEFAULT_SAVE_FILE},
    sequence::seq_to_str,
    snapgene_parse::{export_snapgene, import_snapgene},
    State,
};
use crate::genbank_parse::export_genbank;

fn save_button(
    dialog: &mut FileDialog,
    plasmid_name: &str,
    extension: &str,
    text: &str,
    hover_text: &str,
    ui: &mut Ui,
) {
    if ui.button(text).on_hover_text(hover_text).clicked() {
        let mut save_path = env::current_dir().unwrap();

        let filename = {
            let name = if plasmid_name.is_empty() {
                "a_plasmid".to_string()
            } else {
                plasmid_name.to_lowercase().replace(" ", "_")
            };
            format!("{name}.{extension}")
        };
        save_path.push(Path::new(&filename));

        println!("Save path: {:?}", save_path);

        dialog.save_file();
    }
}

fn load_button(dialog: &mut FileDialog, text: &str, hover_text: &str, ui: &mut Ui) {
    if ui.button(text).on_hover_text(hover_text).clicked() {
        dialog.select_file();
    }
}

/// Ui elements for saving and loading data in various file formats. This includes our own format,
/// FASTA, and (eventually) SnapGene's DNA format.
pub fn save_section(state: &mut State, ui: &mut Ui) {
    if ui
        .button("Save")
        .on_hover_text("Save data. (Ctrl + s)")
        .clicked()
    {
        if let Err(e) = save(DEFAULT_SAVE_FILE, &StateToSave::from_state(state)) {
            println!("Error saving: {e}");
        }
    }

    save_button(
        &mut state.ui.file_dialogs.save,
        &state.plasmid_name,
        "pcad",
        "Save as",
        "Save data in the PlasCAD format, to a file.",
        ui,
    );

    load_button(
        &mut state.ui.file_dialogs.load,
        "Load",
        "Load data in the PlasCAD format, from a file.",
        ui,
    );

    load_button(
        &mut state.ui.file_dialogs.import,
        "Imp FASTA/GenBank/SnapGene",
        "Import a sequence in the FASTA or .dna (SnapGene) formats",
        ui,
    );

    save_button(
        &mut state.ui.file_dialogs.export_fasta,
        &state.plasmid_name,
        "fasta",
        "Exp FASTA",
        "Export the sequence in the FASTA format. This does not include features or primers.",
        ui,
    );

    save_button(
        &mut state.ui.file_dialogs.export_genbank,
        &state.plasmid_name,
        "gb",
        "Exp GenBank",
        "Export the sequence in the GenBank format.",
        ui,
    );

    save_button(
        &mut state.ui.file_dialogs.export_dna,
        &state.plasmid_name,
        "dna",
        "Exp SnapGene",
        "Export data in the .dna (SnapGene) format",
        ui,
    );

    // todo: DRY.
    let ctx = ui.ctx();
    state.ui.file_dialogs.save.update(ctx);

    state.ui.file_dialogs.save.update(ctx);
    state.ui.file_dialogs.load.update(ctx);
    state.ui.file_dialogs.export_fasta.update(ctx);
    state.ui.file_dialogs.export_genbank.update(ctx);
    state.ui.file_dialogs.export_dna.update(ctx);
    state.ui.file_dialogs.import.update(ctx);
    
    let mut sync = false;

    if let Some(path) = state.ui.file_dialogs.import.take_selected() {
        state.ui.file_dialogs.selected = Some(path.to_owned());

        if let Some(extension) = path.extension().and_then(|ext| ext.to_str()) {
            match extension.to_lowercase().as_ref() {
                // Does this work for FASTQ too?
                "fasta" => {
                    if let Ok((seq, id)) = import_fasta(&path) {
                        state.seq = seq;
                        state.plasmid_name = id;
                        sync = true;
                    }
                }
                "dna" => {
                    if let Ok(data) = import_snapgene(&path) {
                        if let Some(v) = data.seq {
                            state.seq = v;
                        }

                        if let Some(v) = data.topology {
                            state.topology = v;
                        }

                        if let Some(v) = data.features {
                            state.features = v;
                        }

                        if let Some(v) = data.primers {
                            state.primer_data = v.into_iter().map(|v| PrimerData::new(v)).collect();

                            state.sync_primer_metrics();
                        }
                        sync = true;
                    }
                }
                "gb" | "gbk" => {
                    // todo: This is repetative with the above.
                    if let Ok(data) = import_genbank(&path) {
                        // if let Some(v) = data.seq {
                        state.seq = data.seq;
                        // }

                        // if let Some(v) = data.topology {
                        state.topology = data.topology;
                        // }

                        // if let Some(v) = data.features {
                        state.features = data.features;
                        // }

                        // todo: No primers?
                        // if let Some(v) = data.primers {
                        //     state.primer_data = v.into_iter().map(|v| PrimerData::new(v)).collect();
                        //
                        //     state.sync_primer_metrics();
                        // }
                        sync = true;
                    }
                }
                _ => {
                    eprintln!("The file to import must be in FASTA, GenBank, or SnapGene format.")
                }
            }
        }
    } else if let Some(path) = state.ui.file_dialogs.save.take_selected() {
        state.ui.file_dialogs.selected = Some(path.to_owned());

        // todo: DRY with above for filename setting.
        let extension = "pcad";
        let filename = {
            let name = if state.plasmid_name.is_empty() {
                "a_plasmid".to_string()
            } else {
                state.plasmid_name.to_lowercase().replace(" ", "_")
            };
            format!("{name}.{extension}")
        };

        if let Err(e) = save(&filename, &StateToSave::from_state(state)) {
            eprintln!("Error saving in PlasCAD format: {:?}", e);
        };
    } else if let Some(path) = state.ui.file_dialogs.load.take_selected() {
        state.ui.file_dialogs.selected = Some(path.to_owned());
        *state = State::load(&path.to_str().unwrap());
        sync = true;
    } else if let Some(path) = state.ui.file_dialogs.export_fasta.take_selected() {
        state.ui.file_dialogs.selected = Some(path.to_owned());

        if let Err(e) = export_fasta(&state.seq, &state.plasmid_name, &path) {
            eprintln!("Error exporting to FASTA: {:?}", e);
        };
    } else if let Some(path) = state.ui.file_dialogs.export_genbank.take_selected() {
        state.ui.file_dialogs.selected = Some(path.to_owned());

        if let Err(e) = export_genbank(
            &state.seq,
            state.topology,
            &state.features,
            &state.primer_data,
            &path,
        ) {
            eprintln!("Error exporting to GenBank: {:?}", e);
        };
    } else if let Some(path) = state.ui.file_dialogs.export_dna.take_selected() {
        state.ui.file_dialogs.selected = Some(path.to_owned());

        if let Err(e) = export_snapgene(
            &state.seq,
            state.topology,
            &state.features,
            &state.primer_data,
            &path,
        ) {
            eprintln!("Error exporting to SnapGene: {:?}", e);
        };
    }

    if sync {
        state.sync_pcr();
        state.sync_primer_metrics();
        state.sync_seq_related(None);
        state.ui.seq_input = seq_to_str(&state.seq);
    }
}
