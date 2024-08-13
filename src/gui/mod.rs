use std::{
    path::{Path, PathBuf},
    str::FromStr,
};

use eframe::{
    egui,
    egui::{
        pos2, Color32, Context, InputState, Key, PointerButton, ScrollArea, TextEdit, Ui,
        ViewportCommand,
    },
    emath::RectTransform,
};
use navigation::Page;
use url::Url;

use crate::{
    file_io::save::{save, StateToSave, DEFAULT_SAVE_FILE},
    gui::{
        primer_qc::primer_details,
        seq_view::{NT_WIDTH_PX, VIEW_AREA_PAD_LEFT, VIEW_AREA_PAD_RIGHT},
    },
    sequence::{Feature, FeatureType, Nucleotide},
    util, Selection, State,
};

mod circle;
mod cloning;
mod feature_overlay;
mod feature_table;
mod metadata;
pub mod navigation;
mod pcr;
mod portions;
mod primer_arrow;
pub mod primer_qc;
mod save;
pub mod seq_view;
pub mod sequence;
// pub for a few consts

pub const WINDOW_WIDTH: f32 = 1300.;
pub const WINDOW_HEIGHT: f32 = 1_000.;

pub const WINDOW_TITLE: &str = "PlasCAD";

pub const ROW_SPACING: f32 = 22.;
pub const COL_SPACING: f32 = 30.;

// Note: This is basically working, but doesn't seem to reflect this scaling factor accurately.
pub const SPLIT_SCREEN_MAX_HEIGHT: f32 = 3.5;

// todo: Move this BLAST stuff A/R.
const NCBI_BLAST_URL: &str = "https://blast.ncbi.nlm.nih.gov/Blast.cgi";

pub fn int_field(val: &mut usize, label: &str, ui: &mut Ui) {
    ui.label(label);
    let mut entry = val.to_string();
    if ui
        .add(TextEdit::singleline(&mut entry).desired_width(40.))
        .changed()
    {
        *val = entry.parse().unwrap_or(0);
    }
}

/// Get a text-representation of the cursor index; a slightly processed version of the raw index.
/// We use this on the sequence and circle views.
pub fn get_cursor_text(cursor_seq_i: Option<usize>, seq_len: usize) -> String {
    match cursor_seq_i {
        Some(p) => {
            if p < seq_len {
                // + 1, as the convention is to use 1-based indexing vice 0.
                (p + 1).to_string()
                // This occurs if the cursor is on the last row, right of the last NT.
            } else {
                String::new()
            }
        }
        None => String::new(),
    }
}

/// Open the web browser to a NCBI-BLAST page, of the sequence of interest.
///
///Example BLAST
/// note: There appears to be a NT limit that will fail most full plastmids when using the GET api.
/// ?PAGE_TYPE=BlastSearch&CMD=Web&LAYOUT=OneWindow&PROGRAM=blastn&MEGABLAST=on&PAGE=Nucleotides&DATABASE=nr
/// &FORMAT_TYPE=HTML&NCBI_GI=on&SHOW_OVERVIEW=on&QUERY=%3Ettt%20%20(43%20..%20905%20%3D%20863%20bp)
/// %0AACTCACTATAGGGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACCGGTATGACTAGTATGGAAGACGCCAAAAACATAAAGAAAGGCCCG
/// GCGCCATTCTATCCGCTGGAAGATGGAACCGCTGGAGAGCAACTGCATAAGGCTATGAAGAGATACGCCCTGGTTCCTGGAACAATTGCTTTTACAGA
/// TGCACATATCGAGGTGGACATCACTTACGCTGAGTACTTCGAAATGTCCGTTCGGTTGGCAGAAGCTATGAAACGATATGGGCTGAATACAAATCACAGA
/// ATCGTCGTATGCAGTGAAAACTCTCTTCAATTCTTTATGCCGGTGTTGGGCGCGTTATTTATCGGAGTTGCAGTTGCGCCCGCGAACGACATTTATAATGA
/// ACGTGAATTGCTCAACAGTATGGGCATTTCGCAGCCTACCGTGGTGTTCGTTTCCAAAAAGGGGTTGCAAAAAATTTTGAACGTGCAAAAAAAGCTCCCAAT
/// CATCCAAAAAATTATTATCATGGATTCTAAAACGGATTACCAGGGATTTCAGTCGATGTACACGTTCGTCACATCTCATCTACCTCCCGGTTTTAATGAATAC
/// GATTTTGTGCCAGAGTCCTTCGATAGGGACAAGACAATTGCACTGATCATGAACTCCTCTGGATCTACTGGTCTGCCTAAAGGTGTCGCTCTGCCTCATAGAACT
/// GCCTGCGTGAGATTCTCGCATGCCAGAGATCCTATTTTTGGCAATCAAATCATTCCGGATACTGCGATTTTAAGTGTTGTTCCATTCCATCACGGTTTTGGAA
/// TGTTTACTACACTCGGATATTTGATATGTGGATTTCGAGTCGTCTTAATGTATAGAT
/// todo: Copy to clipboard  for longer seqs?
fn open_blast(seq: &[Nucleotide]) {
    let params = vec![
        ("PAGE_TYPE", "BlastSearch"),
        ("CMD", "Web"),
        ("LAYOUT", "OneWindow"),
        ("PROGRAM", "blastn"),
        ("MEGABLAST", "on"),
        ("PAGE", "Nucleotides"),
        ("DATABASE", "nr"),
        ("FORMAT_TYPE", "HTML"),
        ("NCBI_GI", "on"),
        ("SHOW_OVERVIEW", "on"),
        ("QUERY", ">ttt  (43 .. 905 = 863 bp)\nACTCACTATAGGGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACCGGTATGACTAGTATGGAAGACGCCAAAAACATAAAGAAAGGCCCGGCGCCATTCTATCCGCTGGAAGATGGAACCGCTGGAGAGCAACTGCATAAGGCTATGAAGAGATACGCCCTGGTTCCTGGAACAATTGCTTTTACAGATGCACATATCGAGGTGGACATCACTTACGCTGAGTACTTCGAAATGTCCGTTCGGTTGGCAGAAGCTATGAAACGATATGGGCTGAATACAAATCACAGAATCGTCGTATGCAGTGAAAACTCTCTTCAATTCTTTATGCCGGTGTTGGGCGCGTTATTTATCGGAGTTGCAGTTGCGCCCGCGAACGACATTTATAATGAACGTGAATTGCTCAACAGTATGGGCATTTCGCAGCCTACCGTGGTGTTCGTTTCCAAAAAGGGGTTGCAAAAAATTTTGAACGTGCAAAAAAAGCTCCCAATCATCCAAAAAATTATTATCATGGATTCTAAAACGGATTACCAGGGATTTCAGTCGATGTACACGTTCGTCACATCTCATCTACCTCCCGGTTTTAATGAATACGATTTTGTGCCAGAGTCCTTCGATAGGGACAAGACAATTGCACTGATCATGAACTCCTCTGGATCTACTGGTCTGCCTAAAGGTGTCGCTCTGCCTCATAGAACTGCCTGCGTGAGATTCTCGCATGCCAGAGATCCTATTTTTGGCAATCAAATCATTCCGGATACTGCGATTTTAAGTGTTGTTCCATTCCATCACGGTTTTGGAATGTTTACTACACTCGGATATTTGATATGTGGATTTCGAGTCGTCTTAATGTATAGAT"),
    ];

    let mut url = Url::parse(NCBI_BLAST_URL).unwrap();
    {
        let mut query_pairs = url.query_pairs_mut();
        for (key, value) in params {
            query_pairs.append_pair(key, value);
        }
    }

    // Open the URL in the default web browser
    if let Err(e) = webbrowser::open(url.as_str()) {
        eprintln!("Failed to open the web browser: {:?}", e);
    }
}

/// Handle an origin change.
fn origin_change(state: &mut State, ui: &mut Ui) {
    if ui.button("Set origin").clicked() {
        state.ui.show_origin_change = !state.ui.show_origin_change;
    }
    if state.ui.show_origin_change {
        ui.horizontal(|ui| {
            ui.label("New origin:");

            let mut entry = state.ui.new_origin.to_string();
            if ui
                .add(TextEdit::singleline(&mut entry).desired_width(40.))
                .changed()
            {
                state.ui.new_origin = entry.parse().unwrap_or(0);
            }

            if ui.button("Set").clicked() {
                util::change_origin(state);
            }

            // ui.add(
            //     // egui::Slider::from_get_set(0.0..=state.generic.seq.len() as f32, |v| {
            //     egui::Slider::new(&mut 0, 0..=state.generic.seq.len(), |v| {
            //         if let Some(v_) = v {
            //
            //             bases_modified.push(basis_i);
            //         }
            //
            //         v
            //         .text("Wt"),
            // );
        });
    }
}

/// Find the index of the smallest feature that contains an index. Index is in our 1-based system.
fn feature_from_index(index: &Option<usize>, features: &[Feature]) -> Option<usize> {
    if let Some(seq_i) = index {
        // If multiple features are in the cursor's region, choose the smallest.
        let mut smallest_feature = 0;
        let mut smallest_feature_size = 99999;

        for (i, feature) in features.iter().enumerate() {
            if feature.feature_type == FeatureType::Source {
                continue; // From GenBank; generally the whole seq.
            }

            if *seq_i > feature.index_range.0 && *seq_i < feature.index_range.1 {
                let feature_size = feature.index_range.1 - feature.index_range.0;
                if feature_size < smallest_feature_size {
                    smallest_feature = i;
                    smallest_feature_size = feature_size;

                    return Some(i);
                }
            }
        }
    }
    None
}

/// Selects a feature, if there is a click in the appropriate canvas; used in both the sequence,
/// and map views.
pub fn select_feature(state: &mut State, from_screen: &RectTransform) {
    if state.ui.click_pending_handle {
        // Don't let clicks out of this canvas remove the selected item.
        if let Some(pos) = state.ui.cursor_pos {
            let pos_rel = from_screen * pos2(pos.0, pos.1);

            if pos_rel.x > 0. && pos_rel.y > 0. {
                let feature_i = feature_from_index(&state.ui.cursor_seq_i, &state.generic.features);

                let mut toggled_off = false;
                if let Selection::Feature(j) = state.ui.selected_item {
                    if j == feature_i.unwrap_or(999) {
                        state.ui.selected_item = Selection::None;
                        toggled_off = true;
                    }
                }

                if !toggled_off {
                    state.ui.selected_item = match feature_i {
                        Some(i) => Selection::Feature(i),
                        None => Selection::None,
                    };
                }
            }
        }

        state.ui.click_pending_handle = false;
    }
}

pub fn set_window_title(path_loaded: &Option<PathBuf>, ui: &mut Ui) {
    // Update the tilebar to reflect the current path loaded or saved.
    if let Some(path) = path_loaded {
        let filename = path
            .file_name()
            .and_then(|name| name.to_str())
            .map(|name_str| name_str.to_string())
            .unwrap();

        let window_title = format!("{WINDOW_TITLE} - {}", filename);
        ui.ctx()
            .send_viewport_cmd(ViewportCommand::Title(window_title));
    }
}

/// Handles keyboard and mouse input not associated with a widget.
/// todo: MOve to a separate module if this becomes complex.
fn handle_input(state: &mut State, ctx: &Context) {
    ctx.input(|ip| {
        if ip.key_pressed(Key::A) && ip.modifiers.ctrl {
            state.generic.primers.push(Default::default());
        }

        if ip.key_pressed(Key::S) && ip.modifiers.ctrl {
            if let Err(e) = save(
                &PathBuf::from(DEFAULT_SAVE_FILE),
                &StateToSave::from_state(state),
            ) {
                println!("Error saving: {e}");
            }
        }

        state.ui.cursor_pos = ip.pointer.hover_pos().map(|pos| (pos.x, pos.y));

        if ip.pointer.button_clicked(PointerButton::Primary) {
            state.ui.click_pending_handle = true;
        }

        // This is a bit awk; borrow errors.
        let mut move_cursor = None;
        // todo: How can we control the rate?
        if let Some(i) = &mut state.ui.text_cursor_i {
            if ip.key_pressed(Key::ArrowLeft) {
                move_cursor = Some(-1);
            }
            if ip.key_pressed(Key::ArrowRight) {
                move_cursor = Some(1);
            }
            if ip.key_pressed(Key::ArrowUp) {
                // todo: DRY with seq_view
                move_cursor = Some(-(state.ui.nt_chars_per_row as isize));
            }
            if ip.key_pressed(Key::ArrowDown) {
                Some(state.ui.nt_chars_per_row);
            }
        }

        if let Some(i) = state.ui.text_cursor_i {
            // Add NTs.
            if ip.key_pressed(Key::A) && !ip.modifiers.ctrl {
                state.insert_nucleotides(&[Nucleotide::A], i);
            }
            if ip.key_pressed(Key::T) {
                state.insert_nucleotides(&[Nucleotide::T], i);
            }
            if ip.key_pressed(Key::C) {
                state.insert_nucleotides(&[Nucleotide::C], i);
            }
            if ip.key_pressed(Key::G) {
                state.insert_nucleotides(&[Nucleotide::G], i);
            }
        }

        if ip.key_pressed(Key::A) && !ip.modifiers.ctrl {
            move_cursor = Some(1);
        }
        if ip.key_pressed(Key::T) {
            move_cursor = Some(1);
        }
        if ip.key_pressed(Key::C) {
            move_cursor = Some(1);
        }
        if ip.key_pressed(Key::G) {
            move_cursor = Some(1);
        }

        if ip.key_pressed(Key::Backspace) {
            // todo: IMpl deleting.
            move_cursor = Some(-1);
        }

        if ip.key_pressed(Key::Delete) {
            move_cursor = Some(-1);
        }

        if let Some(i) = &mut state.ui.text_cursor_i {
            if let Some(amt) = move_cursor {
                // todo: Impl limits
                *i += amt as usize;
            }

            // if decrement_cursor {
            //     if *i > 0 {
            //         *i -= 1;
            //     }
            // } else if increment_cursor {
            //     if *i + 1 < state.generic.seq.len() {
            //         *i += 1;
            //     }
            // }
        }
    });
}

pub fn draw(state: &mut State, ctx: &Context) {
    handle_input(state, ctx);

    egui::CentralPanel::default().show(ctx, |ui| {
        // todo: This section DRY with seq viewx.

        let mut visuals = ctx.style().visuals.clone();
        // visuals.override_text_color = Some(Color32::from_rgb(255, 0, 0));
        visuals.override_text_color = Some(Color32::LIGHT_GRAY);
        ctx.set_visuals(visuals);

        ui.horizontal(|ui| {
            navigation::page_selector(state, ui);

            ui.add_space(COL_SPACING);

            ui.label("Name: ");
            ui.add(
                TextEdit::singleline(&mut state.generic.metadata.plasmid_name).desired_width(280.),
            );
        });

        ui.add_space(ROW_SPACING / 2.);

        ui.horizontal(|ui| {
            save::save_section(state, ui);

            ui.add_space(COL_SPACING);

            // todo: YOu will need a better organization method.
            if ui.button("BLAST").clicked() {
                open_blast(&state.generic.seq); // todo: Seq A/R
            }

            origin_change(state, ui);
        });

        ui.add_space(ROW_SPACING / 2.);

        ScrollArea::vertical().show(ui, |ui| match state.ui.page {
            Page::Sequence => sequence::seq_page(state, ui),
            Page::Map => circle::circle_page(state, ui),
            Page::Features => feature_table::features_page(state, ui),
            Page::Primers => primer_details(state, ui),
            Page::Cloning => {
                cloning::seq_editor_slic(state, ui);
            }
            Page::Pcr => pcr::pcr_page(state, ui),
            Page::Metadata => metadata::metadata_page(&mut state.generic.metadata, ui),
            _ => (),
            // Page::Portions => portions::portions_page(state, ui),
        });
    });
}
