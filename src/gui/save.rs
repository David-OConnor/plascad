//! GUI code for saving and loading. Calls business logic in `file_io/save.rs`.

use std::{env, path::Path};

use eframe::egui::Ui;
use egui_file_dialog::FileDialog;

use crate::{
    file_io::{
        genbank::export_genbank,
        save,
        save::{export_fasta, load_import, StateToSave},
        snapgene::export_snapgene,
    },
    gui::{navigation::Tab, set_window_title},
    state::State,
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
        dialog.pick_file();
    }
}

/// Ui elements for saving and loading data in various file formats. This includes our own format,
/// FASTA, and (eventually) SnapGene's DNA format.
pub fn save_section(state: &mut State, ui: &mut Ui) {
    let button_text = if state.tabs_open[state.active].path.is_some() {
        "Save"
    } else {
        "Quicksave"
    };
    if ui
        .button(button_text)
        .on_hover_text("Save data. (Ctrl + S)")
        .clicked()
    {
        save::save_current_file(state);
    }

    save_button(
        &mut state.ui.file_dialogs.save,
        &state.generic[state.active].metadata.plasmid_name,
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
        &state.generic[state.active].metadata.plasmid_name,
        "fasta",
        "Exp FASTA",
        "Export the sequence in the FASTA format. This does not include features or primers.",
        ui,
    );

    save_button(
        &mut state.ui.file_dialogs.export_genbank,
        &state.generic[state.active].metadata.plasmid_name,
        "gbk",
        "Exp GenBank",
        "Export data in the GenBank format.",
        ui,
    );

    save_button(
        &mut state.ui.file_dialogs.export_dna,
        &state.generic[state.active].metadata.plasmid_name,
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

    if let Some(path) = state.ui.file_dialogs.load.take_picked() {
        sync = true;
        if let Some(loaded) = load_import(&path) {
            state.load(&loaded);
        }
    } else if let Some(path) = state.ui.file_dialogs.save.take_picked() {
        match StateToSave::from_state(state, state.active).save_to_file(&path) {
            Ok(_) => {
                state.tabs_open[state.active] = Tab {
                    path: Some(path.to_owned()),
                    ab1: false,
                };
                set_window_title(&state.tabs_open[state.active], ui);
                state.save_prefs(); // Save opened tabs.
            }
            Err(e) => eprintln!("Error saving in PlasCAD format: {:?}", e),
        };
    } else if let Some(path) = state.ui.file_dialogs.export_fasta.take_picked() {
        match export_fasta(
            state.get_seq(),
            &state.generic[state.active].metadata.plasmid_name,
            &path,
        ) {
            Ok(_) => {
                state.tabs_open[state.active] = Tab {
                    path: Some(path.to_owned()),
                    ab1: false,
                };
                set_window_title(&state.tabs_open[state.active], ui);
                state.save_prefs(); // Save opened tabs.
            }
            Err(e) => eprintln!("Error exporting to FASTA: {:?}", e),
        }
    } else if let Some(path) = state.ui.file_dialogs.export_genbank.take_picked() {
        let mut primer_matches = Vec::new();
        for primer in &state.generic[state.active].primers {
            for prim_match in &primer.volatile.matches {
                primer_matches.push((prim_match.clone(), primer.name.clone()));
            }
        }

        match export_genbank(&state.generic[state.active], &primer_matches, &path) {
            Ok(_) => {
                state.tabs_open[state.active] = Tab {
                    path: Some(path.to_owned()),
                    ab1: false,
                };
                set_window_title(&state.tabs_open[state.active], ui);
                state.save_prefs(); // Save opened tabs.
            }
            Err(e) => eprintln!("Error exporting to GenBank: {:?}", e),
        }
    } else if let Some(path) = state.ui.file_dialogs.export_dna.take_picked() {
        match export_snapgene(&state.generic[state.active], &path) {
            Ok(_) => {
                state.tabs_open[state.active] = Tab {
                    path: Some(path.to_owned()),
                    ab1: false,
                };
                set_window_title(&state.tabs_open[state.active], ui);
                state.save_prefs(); // Save opened tabs.
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

        set_window_title(&state.tabs_open[state.active], ui);
    }
}
