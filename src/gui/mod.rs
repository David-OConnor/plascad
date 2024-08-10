use std::{path::PathBuf, str::FromStr};

use eframe::{
    egui,
    egui::{Color32, Context, Key, ScrollArea, TextEdit, Ui},
};
use eframe::egui::PointerButton;
use navigation::Page;
use url::Url;

use crate::{
    file_io::save::{save, StateToSave, DEFAULT_SAVE_FILE},
    gui::primer_qc::primer_details,
    sequence::{Feature, FeatureType, Nucleotide},
    util, State,
};

mod circle;
mod feature_overlay;
mod features;
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
                TextEdit::singleline(&mut state.generic.metadata.plasmid_name).desired_width(200.),
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

        ui.add_space(ROW_SPACING);

        ScrollArea::vertical().show(ui, |ui| match state.ui.page {
            Page::Sequence => sequence::seq_page(state, ui),
            Page::Map => circle::circle_page(state, ui),
            Page::Features => features::features_page(state, ui),
            Page::Primers => primer_details(state, ui),
            Page::Pcr => pcr::pcr_page(state, ui),
            Page::Metadata => metadata::metadata_page(&mut state.generic.metadata, ui),
            _ => (),
            // Page::Portions => portions::portions_page(state, ui),
        });
    });
}
