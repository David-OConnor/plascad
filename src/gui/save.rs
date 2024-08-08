//! GUI code for saving and loading

use std::{
    env,
    ops::Range,
    path::{Path, PathBuf},
};

use eframe::egui::Ui;
// use egui_file::FileDialog;
use egui_file_dialog::FileDialog;

use crate::{
    file_io::{
        genbank::{export_genbank, import_genbank},
        save::{export_fasta, import_fasta, save, StateToSave, DEFAULT_SAVE_FILE},
        snapgene::{export_snapgene, import_snapgene},
    },
    primer::{PrimerData, PrimerDirection},
    sequence::seq_to_str,
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
                plasmid_name.to_lowercase().replace(" ", "_")
            };
            format!("{name}.{extension}")
        };
        save_path.push(Path::new(&filename));

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
        if let Err(e) = save(
            &PathBuf::from(DEFAULT_SAVE_FILE),
            &StateToSave::from_state(state),
        ) {
            println!("Error saving: {e}");
        }
    }

    save_button(
        &mut state.ui.file_dialogs.save,
        &state.generic.metadata.plasmid_name,
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
        &state.generic.metadata.plasmid_name,
        "fasta",
        "Exp FASTA",
        "Export the sequence in the FASTA format. This does not include features or primers.",
        ui,
    );

    save_button(
        &mut state.ui.file_dialogs.export_genbank,
        &state.generic.metadata.plasmid_name,
        "gb",
        "Exp GenBank",
        "Export the sequence in the GenBank format.",
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
    state.ui.file_dialogs.import.update(ctx);

    let mut sync = false;

    if let Some(path) = state.ui.file_dialogs.import.take_selected() {
        state.ui.file_dialogs.selected = Some(path.to_owned());

        if let Some(extension) = path.extension().and_then(|ext| ext.to_str()) {
            match extension.to_lowercase().as_ref() {
                // Does this work for FASTQ too?
                "fasta" => {
                    if let Ok((seq, id, description)) = import_fasta(&path) {
                        state.generic.seq = seq;
                        state.generic.metadata.plasmid_name = id;
                        state.generic.metadata.comments = vec![description];
                        sync = true;
                    }
                }
                "dna" => {
                    if let Ok(data) = import_snapgene(&path) {
                        state.generic = data;
                        sync = true;
                    }
                }
                "gb" | "gbk" => {
                    if let Ok(data) = import_genbank(&path) {
                        state.generic = data;
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
        if let Err(e) = save(&path, &StateToSave::from_state(state)) {
            eprintln!("Error saving in PlasCAD format: {:?}", e);
        };
    } else if let Some(path) = state.ui.file_dialogs.load.take_selected() {
        state.ui.file_dialogs.selected = Some(path.to_owned());
        *state = State::load(&path.to_str().unwrap());
        sync = true;
    } else if let Some(path) = state.ui.file_dialogs.export_fasta.take_selected() {
        state.ui.file_dialogs.selected = Some(path.to_owned());

        if let Err(e) = export_fasta(
            &state.generic.seq,
            &state.generic.metadata.plasmid_name,
            &path,
        ) {
            eprintln!("Error exporting to FASTA: {:?}", e);
        };
    } else if let Some(path) = state.ui.file_dialogs.export_genbank.take_selected() {
        state.ui.file_dialogs.selected = Some(path.to_owned());

        let mut primer_matches = Vec::new();
        for primer in &state.generic.primers {
            for (dir, range) in &primer.volatile.matches_seq {
                primer_matches.push((*dir, range.clone(), primer.name.clone()));
            }
        }

        if let Err(e) = export_genbank(&state.generic, &primer_matches, &path) {
            eprintln!("Error exporting to GenBank: {:?}", e);
        };
    } else if let Some(path) = state.ui.file_dialogs.export_dna.take_selected() {
        state.ui.file_dialogs.selected = Some(path.to_owned());

        if let Err(e) = export_snapgene(&state.generic, &path) {
            eprintln!("Error exporting to SnapGene: {:?}", e);
        };
    }

    if sync {
        state.sync_pcr();
        state.sync_primer_metrics();
        state.sync_seq_related(None);
        state.ui.seq_input = seq_to_str(&state.generic.seq);
    }
}
