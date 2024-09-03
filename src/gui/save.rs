//! GUI code for saving and loading

use std::{
    env,
    path::{Path, PathBuf},
    str::FromStr,
};

use eframe::egui::Ui;
use egui_file_dialog::FileDialog;

use crate::{
    feature_db_load::find_features,
    file_io::{
        genbank::{export_genbank, import_genbank},
        save::{
            export_fasta, import_fasta, save, StateToSave, DEFAULT_PREFS_FILE, DEFAULT_SAVE_FILE,
        },
        snapgene::{export_snapgene, import_snapgene},
    },
    gui::set_window_title,
    State,
};

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
                plasmid_name.to_lowercase().replace(' ', "_")
            };
            format!("{name}.{extension}")
        };
        save_path.push(Path::new(&filename));

        dialog.config_mut().default_file_name = filename.to_string();
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
        .button("New")
        .on_hover_text("Create and open a new file. (Ctrl + N)")
        .clicked()
    {
        state.reset();
        set_window_title(&state.path_loaded, ui);
    }

    let button_text = if state.path_loaded.is_some() {
        "Save"
    } else {
        "Quicksave"
    };
    if ui
        .button(button_text)
        .on_hover_text("Save data. (Ctrl + S)")
        .clicked()
    {
        save_current_file(state);
    }

    save_button(
        &mut state.ui.file_dialogs.save,
        &state.generic.metadata.plasmid_name,
        "pcad",
        "Save as",
        "Save data in the PlasCAD format. (Ctrl + Shift + S)",
        ui,
    );

    load_button(
        &mut state.ui.file_dialogs.load,
        "Load/Import",
        "Load data in the PlasCAD, FASTA, GenBank, or .dna (SnapGene) formats (Ctrl + O)",
        ui,
    );

    save_button(
        &mut state.ui.file_dialogs.export_fasta,
        &state.generic.metadata.plasmid_name,
        "fasta",
        "Exp FASTA",
        "Export the sequence in the FASTA format. This does not include features or primers.",
        ui,
    );

    save_button(
        &mut state.ui.file_dialogs.export_genbank,
        &state.generic.metadata.plasmid_name,
        "gbk",
        "Exp GenBank",
        "Export data in the GenBank format.",
        ui,
    );

    save_button(
        &mut state.ui.file_dialogs.export_dna,
        &state.generic.metadata.plasmid_name,
        "dna",
        "Exp SnapGene",
        "Export data in the .dna (SnapGene) format",
        ui,
    );

    // todo: DRY.
    let ctx = ui.ctx();

    state.ui.file_dialogs.save.update(ctx);
    state.ui.file_dialogs.load.update(ctx);
    state.ui.file_dialogs.export_fasta.update(ctx);
    state.ui.file_dialogs.export_genbank.update(ctx);
    state.ui.file_dialogs.export_dna.update(ctx);

    let mut sync = false;

    if let Some(path) = state.ui.file_dialogs.load.take_selected() {
        sync = true;
        load_import(state, &path);
    } else if let Some(path) = state.ui.file_dialogs.save.take_selected() {
        match save(&path, &StateToSave::from_state(state)) {
            Ok(_) => {
                state.path_loaded = Some(path.to_owned());
                state.save_prefs(); // to sync last opened.

                set_window_title(&state.path_loaded, ui);
            }
            Err(e) => eprintln!("Error saving in PlasCAD format: {:?}", e),
        };
    } else if let Some(path) = state.ui.file_dialogs.export_fasta.take_selected() {
        match export_fasta(
            &state.generic.seq,
            &state.generic.metadata.plasmid_name,
            &path,
        ) {
            Ok(_) => {
                state.path_loaded = Some(path.to_owned());
                set_window_title(&state.path_loaded, ui);
            }
            Err(e) => eprintln!("Error exporting to FASTA: {:?}", e),
        }
    } else if let Some(path) = state.ui.file_dialogs.export_genbank.take_selected() {
        let mut primer_matches = Vec::new();
        for primer in &state.generic.primers {
            for prim_match in &primer.volatile.matches {
                primer_matches.push((prim_match.clone(), primer.name.clone()));
            }
        }

        match export_genbank(&state.generic, &primer_matches, &path) {
            Ok(_) => {
                state.path_loaded = Some(path.to_owned());
                set_window_title(&state.path_loaded, ui);
            }
            Err(e) => eprintln!("Error exporting to GenBank: {:?}", e),
        }
    } else if let Some(path) = state.ui.file_dialogs.export_dna.take_selected() {
        match export_snapgene(&state.generic, &path) {
            Ok(_) => {
                state.path_loaded = Some(path.to_owned());
                set_window_title(&state.path_loaded, ui);
            }
            Err(e) => eprintln!("Error exporting to SnapGene: {:?}", e),
        };
    }

    if sync {
        state.sync_pcr();
        state.sync_primer_metrics();
        state.sync_seq_related(None);
        state.sync_portions();
        state.reset_selections();

        set_window_title(&state.path_loaded, ui);
    }
}

/// Load state from a file of various formats.
pub fn load_import(state: &mut State, path: &Path) {
    if let Some(extension) = path.extension().and_then(|ext| ext.to_str()) {
        match extension.to_lowercase().as_ref() {
            "pcad" => {
                *state = State::load(&path, &PathBuf::from_str(DEFAULT_PREFS_FILE).unwrap());
                state.path_loaded = Some(path.to_owned());
                state.save_prefs(); // to sync last opened.
            }
            // Does this work for FASTQ too?
            "fasta" => {
                if let Ok((seq, id, description)) = import_fasta(&path) {
                    state.generic.seq = seq;
                    state.generic.metadata.plasmid_name = id;
                    state.generic.metadata.comments = vec![description];
                    // FASTA is seq-only data, so don't attempt to save over it.
                    state.path_loaded = None;
                    state.save_prefs();

                    // Automatically annotate FASTA files.
                    state.generic.features = find_features(&state.generic.seq);
                }
            }
            "dna" => {
                if let Ok(data) = import_snapgene(&path) {
                    state.generic = data;
                    // We do not mark the path as opened if using SnapGene, since we currently can not
                    // fully understand the format, nor make a native file SnapGene can open.
                    state.path_loaded = None;
                }
            }
            "gb" | "gbk" => {
                if let Ok(data) = import_genbank(&path) {
                    state.generic = data;
                    state.path_loaded = Some(path.to_owned());
                    state.save_prefs();
                }
            }
            _ => {
                eprintln!(
                    "The file to import must be in PlasCAD, FASTA, GenBank, or SnapGene format."
                )
            }
        }
    }
}

/// Save the current file ("save" vice "save as") if there is one; if not, quicksave to an anonymous file.
pub fn save_current_file(state: &State) {
    match &state.path_loaded {
        Some(path) => {
            if let Some(extension) = path.extension().and_then(|ext| ext.to_str()) {
                match extension.to_lowercase().as_ref() {
                    "pcad" => {
                        if let Err(e) = save(path, &StateToSave::from_state(state)) {
                            eprintln!("Error saving in PlasCAD format: {:?}", e);
                        };
                    }
                    // Does this work for FASTQ too?
                    "fasta" => {
                        if let Err(e) = export_fasta(
                            &state.generic.seq,
                            &state.generic.metadata.plasmid_name,
                            path,
                        ) {
                            eprintln!("Error exporting to FASTA: {:?}", e);
                        };
                    }
                    "dna" => {
                        if let Err(e) = export_snapgene(&state.generic, path) {
                            eprintln!("Error exporting to SnapGene: {:?}", e);
                        };
                    }
                    "gb" | "gbk" => {
                        let mut primer_matches = Vec::new();
                        for primer in &state.generic.primers {
                            for prim_match in &primer.volatile.matches {
                                primer_matches.push((prim_match.clone(), primer.name.clone()));
                            }
                        }

                        if let Err(e) = export_genbank(&state.generic, &primer_matches, path) {
                            eprintln!("Error exporting to GenBank: {:?}", e);
                        };
                    }
                    _ => {
                        eprintln!("Unexpected file format loading.")
                    }
                }
            }
        }
        None => {
            // Quicksave.
            if let Err(e) = save(
                &PathBuf::from(DEFAULT_SAVE_FILE),
                &StateToSave::from_state(state),
            ) {
                eprintln!("Error quicksaving: {e}");
            }
        }
    }
    // todo: You will likely more this to an automatic one.
    state.save_prefs()
}
