//! GUI entry point
//!
//! [Useful emojis supported by EGUI](https://github.com/emilk/egui/blob/9a1e358a144b5d2af9d03a80257c34883f57cf0b/crates/egui/src/lib.rs#L557-L575)
//! ∞⊗⎗⎘⎙⏏⏴⏵⏶⏷
//! ⏩⏪⏭⏮⏸⏹⏺■▶📾🔀🔁🔃
//! ☀☁★☆☐☑☜☝☞☟⛃⛶✔
//! ↺↻⟲⟳⬅➡⬆⬇⬈⬉⬊⬋⬌⬍⮨⮩⮪⮫
//! ♡
//! 📅📆
//! 📈📉📊
//! 📋📌📎📤📥🔆
//! 🔈🔉🔊🔍🔎🔗🔘
//! 🕓🖧🖩🖮🖱🖴🖵🖼🗀🗁🗋🗐🗑🗙🚫❓
//!
//! Also maybe: http://jslegers.github.io/emoji-icon-font/
//! https://fonts.google.com/noto/specimen/Noto+Emoji
//!
//! Note: business logic should be kept out of this module (and its sub-modules), when practical.

use std::path::PathBuf;

use eframe::{
    egui,
    egui::{pos2, Color32, Context, RichText, TextEdit, ThemePreference, Ui, ViewportCommand},
    emath::RectTransform,
};
use navigation::Page;

use crate::{
    external_websites,
    feature_db_load::find_features,
    gui::{
        input::handle_input, navigation::Tab, primer_table::primer_details, theme::COLOR_ACTION,
    },
    misc_types::{Feature, FeatureType},
    primer::Primer,
    state::State,
    util,
    util::{get_window_title, merge_feature_sets},
    Selection,
};

mod ab1;
mod alignment;
mod circle;
mod cloning;
mod feature_table;
mod input;
mod ligation;
mod lin_maps;
mod metadata;
pub mod navigation;
mod pcr;
mod portions;
pub mod primer_table;
mod protein;
pub mod save;
pub mod sequence;
mod theme;

pub const WINDOW_WIDTH: f32 = 1300.;
pub const WINDOW_HEIGHT: f32 = 1_000.;

pub const WINDOW_TITLE: &str = "PlasCAD";

pub const ROW_SPACING: f32 = 22.;
pub const COL_SPACING: f32 = 30.;

// Note: This is basically working, but doesn't seem to reflect this scaling factor accurately.
pub const SPLIT_SCREEN_MAX_HEIGHT: f32 = 3.5;

pub const PRIMER_FWD_COLOR: Color32 = Color32::from_rgb(255, 0, 255);
pub const PRIMER_REV_COLOR: Color32 = Color32::LIGHT_YELLOW;

pub const BACKGROUND_COLOR: Color32 = Color32::from_rgb(10, 20, 10);

pub const COLOR_SEQ: Color32 = Color32::LIGHT_BLUE;
// (0xAD, 0xD8, 0xE6)
pub const COLOR_SEQ_DIMMED: Color32 = Color32::from_rgb(140, 160, 165); // Eg dim when there are search results
pub const COLOR_RE: Color32 = Color32::LIGHT_RED;
pub const COLOR_RE_HIGHLIGHTED: Color32 = Color32::YELLOW;

// If using a dedicated canvas for a linear map.
pub const LINEAR_MAP_HEIGHT: f32 = 60.;

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

/// Get a text-representation of the cursor index (Mouse or text); a slightly processed version of the raw index.
/// We use this on the sequence and circle views.
pub fn get_cursor_text(cursor_seq_i: Option<usize>, seq_len: usize) -> String {
    match cursor_seq_i {
        Some(p) => {
            if p <= seq_len {
                p.to_string()
                // This occurs if the cursor is on the last row, right of the last NT.
            } else {
                String::new()
            }
        }
        None => String::new(),
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

            if ui
                .button(RichText::new("Set").color(COLOR_ACTION))
                .clicked()
            {
                util::change_origin(state);
            }

            // ui.add(
            //     // egui::Slider::from_get_set(0.0..=state.get_seq().len() as f32, |v| {
            //     egui::Slider::new(&mut 0, 0..=state.get_seq().len(), |v| {
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
        let mut matched = false;
        let mut smallest_feature = 0;
        let mut smallest_feature_size = 99999;

        for (i, feature) in features.iter().enumerate() {
            if feature.feature_type == FeatureType::Source {
                continue; // From GenBank; generally the whole seq.
            }

            if *seq_i > feature.range.start && *seq_i < feature.range.end {
                matched = true;
                let feature_size = feature.range.end - feature.range.start - 1;
                if feature_size < smallest_feature_size {
                    smallest_feature = i;
                    smallest_feature_size = feature_size;
                }
            }
        }
        if matched {
            return Some(smallest_feature);
        }
    }
    None
}

/// todo: DRY with `feature_from_index`. Combine.
fn primer_from_index(index: &Option<usize>, primers: &[Primer]) -> Option<usize> {
    if let Some(seq_i) = index {
        // If multiple features are in the cursor's region, choose the smallest.
        let mut matched = false;
        let mut smallest_feature = 0;
        let mut smallest_feature_size = 99999;

        for (i, primer) in primers.iter().enumerate() {
            for p_match in &primer.volatile.matches {
                if *seq_i > p_match.range.start && *seq_i < p_match.range.end {
                    matched = true;
                    let feature_size = p_match.range.end - p_match.range.start - 1;
                    if feature_size < smallest_feature_size {
                        smallest_feature = i;
                        smallest_feature_size = feature_size;
                    }
                }
            }
        }
        if matched {
            return Some(smallest_feature);
        }
    }
    None
}

/// Selects a primer or feature, if there is a click in the appropriate canvas; used in both the sequence,
/// and map views. Currently, if there is a primer and fetaure overlayed, it selects the primer. (This isn't ideal,
/// but acceptable for now.)
pub fn select_feature(state: &mut State, from_screen: &RectTransform) {
    let click_handle = match &mut state.ui.page {
        Page::Sequence => &mut state.ui.dblclick_pending_handle,
        Page::Map => &mut state.ui.click_pending_handle,
        _ => &mut false,
    };

    if *click_handle {
        // Don't let clicks out of this canvas remove the selected item.
        if let Some(pos) = state.ui.cursor_pos {
            let pos_rel = from_screen * pos2(pos.0, pos.1);

            if pos_rel.x > 0. && pos_rel.y > 0. {
                let feature_i = feature_from_index(
                    &state.ui.cursor_seq_i,
                    &state.generic[state.active].features,
                );
                let primer_i =
                    primer_from_index(&state.ui.cursor_seq_i, &state.generic[state.active].primers);

                let mut toggled_off = false;
                if let Selection::Feature(j) = state.ui.selected_item {
                    if primer_i.is_none() {
                        // If primer_i is some, we select it vice selecting None.
                        if j == feature_i.unwrap_or(999) {
                            state.ui.selected_item = Selection::None;
                            toggled_off = true;
                        }
                    }
                }

                // todo: DRY
                if let Selection::Primer(j) = state.ui.selected_item {
                    if j == primer_i.unwrap_or(999) {
                        state.ui.selected_item = Selection::None;
                        toggled_off = true;
                    }
                }

                if !toggled_off {
                    state.ui.selected_item = if let Some(v) = primer_i {
                        Selection::Primer(v)
                    } else {
                        match feature_i {
                            Some(i) => Selection::Feature(i),
                            None => Selection::None,
                        }
                    }
                }
            }
        }

        *click_handle = false;
    }
}

/// Update the tilebar to reflect the current path loaded or saved.
pub fn set_window_title(tab: &Tab, ui: &mut Ui) {
    let title = match &tab.path {
        Some(path) => get_window_title(&path),
        None => WINDOW_TITLE.to_owned(),
    };

    ui.ctx().send_viewport_cmd(ViewportCommand::Title(title));
}

pub fn draw(state: &mut State, ctx: &Context) {
    ctx.options_mut(|o| o.theme_preference = ThemePreference::Dark);

    egui::CentralPanel::default().show(ctx, |ui| {
        handle_input(state, ui);

        // todo: This section DRY with seq viewx.

        let mut visuals = ctx.style().visuals.clone();
        // visuals.override_text_color = Some(Color32::from_rgb(255, 0, 0));
        visuals.override_text_color = Some(Color32::LIGHT_GRAY);
        ctx.set_visuals(visuals);

        navigation::tab_selector(state, ui);
        ui.add_space(ROW_SPACING / 2.);

        let ab1_mode = state.tabs_open[state.active].ab1;

        if ab1_mode {
            // Ensure that if an AB1 file is selected, we show the AB1 page and no tabs.
            // Otherwise, ensure we do not show the AB1 page.
            state.ui.page = Page::Ab1;
        } else {
            if state.ui.page == Page::Ab1 {
                state.ui.page = Page::Map;
            }

            ui.horizontal(|ui| {
                navigation::page_selector(state, ui);

                ui.add_space(COL_SPACING / 2.);

                ui.label("Name:");
                ui.add(
                    TextEdit::singleline(&mut state.generic[state.active].metadata.plasmid_name)
                        .desired_width(260.),
                );

                ui.label(format!("{} bp", state.get_seq().len()));
            });
        }

        ui.add_space(ROW_SPACING / 2.);

        ui.horizontal(|ui| {
            save::save_section(state, ui);

            ui.add_space(COL_SPACING);

            origin_change(state, ui);

            if ui.button("Annotate").clicked() {
                // Don't add duplicates.
                let features = find_features(&state.get_seq());
                merge_feature_sets(&mut state.generic[state.active].features, &features)
            }

            // todo: Kludge.
            if ui.button("Sync RE sites").clicked() {
                state.sync_re_sites();
            }

            ui.add_space(COL_SPACING);
            ui.label("Edit lock:");
            let (lock_text, lock_color) = if state.ui.seq_edit_lock {
                ("🔒", Color32::LIGHT_BLUE)
            } else {
                ("🔓", Color32::from_rgb(255, 210, 140))
            };
            // todo: IDeally bigger font size, but without making the whole line take up more vertical space.
            if ui
                .button(RichText::new(lock_text).color(lock_color))
                .on_hover_text("Prevent edits to the sequence")
                .clicked()
            {
                state.ui.seq_edit_lock = !state.ui.seq_edit_lock;
            }

            // todo: YOu will need a better organization method.
            if state.ui.text_selection.is_some() || state.ui.selected_item != Selection::None {
                let text = if state.ui.text_selection.is_some() {
                    "BLAST selection"
                } else {
                    match state.ui.selected_item {
                        Selection::Feature(_) => "BLAST feature",
                        Selection::Primer(_) => "BLAST primer",
                        Selection::None => unreachable!(),
                    }
                };

                if ui.button(RichText::new(text).color(COLOR_ACTION)).clicked() {
                    external_websites::blast(state);
                }
            }

            ui.add_space(COL_SPACING);

            let mut selection_avail = false;
            if state.ui.text_selection.is_some() {
                selection_avail = true;
            }
            match state.ui.selected_item {
                Selection::None => (),
                _ => {
                    selection_avail = true;
                }
            }

            if selection_avail {
                ui.add_space(COL_SPACING);
                if ui
                    .button("🗐")
                    .on_hover_text("Copy the selected selection, feature or primer. (Ctrl + C)")
                    .clicked()
                {
                    state.copy_seq()
                }
            }
        });

        ui.add_space(ROW_SPACING / 2.);

        // ScrollArea::vertical().show(ui, |ui| match state.ui.page {
        match state.ui.page {
            Page::Sequence => sequence::seq_page(state, ui),
            Page::Map => circle::circle_page(state, ui),
            Page::Features => feature_table::features_page(state, ui),
            Page::Primers => primer_details(state, ui),
            Page::Pcr => pcr::pcr_page(state, ui),
            Page::Alignment => alignment::alignment_page(state, ui),
            Page::Cloning => cloning::cloning_page(state, ui),
            Page::Proteins => protein::protein_page(state, ui),
            Page::Ligation => ligation::ligation_page(state, ui),
            Page::Metadata => {
                metadata::metadata_page(&mut state.generic[state.active].metadata, ui)
            }
            Page::Portions => portions::portions_page(&mut state.portions[state.active], ui),
            Page::Ab1 => {
                if state.ab1_data.len() > 0 {
                    // todo: Fix index, and this hack in general.
                    ab1::ab1_page(&state.ab1_data[0], ui)
                }
            }
        }
    });
}

pub fn select_color_text(text: &str, selected: bool) -> RichText {
    let color = if selected {
        Color32::GREEN
    } else {
        Color32::WHITE
    };

    RichText::new(text).color(color)
}
