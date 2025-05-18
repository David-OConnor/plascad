//! GUI code for the features editor and related.

use eframe::egui::{
    Color32, ComboBox, CursorIcon, Frame, RichText, ScrollArea, Stroke, TextEdit, Ui,
};

use crate::{
    Color, Selection,
    gui::{COL_SPACING, ROW_SPACING, int_field, theme::COLOR_ACTION},
    misc_types::{
        Feature,
        FeatureDirection::{self, Forward, Reverse},
        FeatureType,
    },
    state::State,
    util::RangeIncl,
};

const LABEL_EDIT_WIDTH: f32 = 140.;

/// A color selector for use with feature addition and editing.
fn color_picker(val: &mut Option<Color>, feature_color: Color, ui: &mut Ui) {
    let mut color_override = val.is_some();
    if ui.checkbox(&mut color_override, "").changed() {
        if color_override {
            // Default to the feature color when checking.
            *val = Some(feature_color);
        } else {
            *val = None;
            return;
        }
    }

    // Only show the color picker if choosing to override.
    if val.is_none() {
        return;
    }

    let (r, g, b) = val.unwrap();
    let mut color = Color32::from_rgb(r, g, b);
    if ui.color_edit_button_srgba(&mut color).changed() {
        *val = Some((color.r(), color.g(), color.b()));
    }
}

/// A selector for use with feature addition and editing.
/// todo: Generic selector creator?
pub fn direction_picker(val: &mut FeatureDirection, id: usize, ui: &mut Ui) {
    ComboBox::from_id_salt(id)
        .width(74.)
        .selected_text(val.to_string())
        .show_ui(ui, |ui| {
            for dir in [FeatureDirection::None, Forward, Reverse] {
                ui.selectable_value(val, dir, dir.to_string());
            }
        });
}

/// A selector for use with feature addition and editing.
/// todo: Generic selector creator?
fn feature_type_picker(val: &mut FeatureType, id: usize, ui: &mut Ui) {
    ComboBox::from_id_salt(id)
        .width(140.)
        .selected_text(val.to_string())
        .show_ui(ui, |ui| {
            for feature_type in [
                FeatureType::Generic,
                FeatureType::CodingRegion,
                FeatureType::Ori,
                // FeatureType::RnaPolyBindSite,
                FeatureType::RibosomeBindSite,
                FeatureType::AntibioticResistance,
                FeatureType::LongTerminalRepeat,
                FeatureType::Exon,
                FeatureType::Transcript,
                // todo: Source?
            ] {
                ui.selectable_value(val, feature_type, feature_type.to_string());
            }
        });
}

pub fn feature_table(state: &mut State, ui: &mut Ui) {
    feature_add_disp(state, ui);
    ui.add_space(ROW_SPACING);

    let mut removed = None;
    for (i, feature) in state.generic[state.active].features.iter_mut().enumerate() {
        let mut border_width = 0.;
        if let Selection::Feature(j) = state.ui.selected_item {
            if i == j {
                border_width = 1.;
            }
        }

        Frame::none()
            .stroke(Stroke::new(border_width, Color32::LIGHT_RED))
            .inner_margin(border_width)
            .show(ui, |ui| {
                if ui
                    .heading(RichText::new(&feature.label).color(COLOR_ACTION))
                    .on_hover_cursor(CursorIcon::PointingHand)
                    .clicked()
                {
                    state.ui.selected_item = Selection::Feature(i);
                }

                ui.horizontal(|ui| {
                    int_field(&mut feature.range.start, "Start:", ui);
                    int_field(&mut feature.range.end, "End:", ui);

                    ui.label("Label:");
                    ui.add(
                        TextEdit::singleline(&mut feature.label).desired_width(LABEL_EDIT_WIDTH),
                    );

                    ui.label("Type:");
                    feature_type_picker(&mut feature.feature_type, 100 + i, ui);

                    ui.label("Dir:");
                    direction_picker(&mut feature.direction, 300 + i, ui);

                    ui.label("Custom color:");
                    color_picker(
                        &mut feature.color_override,
                        feature.feature_type.color(),
                        ui,
                    );

                    if ui.button("Add note").clicked() {
                        feature.notes.push((String::new(), String::new()));
                    }

                    // todo: This section repetative with primers.
                    let mut selected = false;
                    if let Selection::Feature(sel_i) = state.ui.selected_item {
                        if sel_i == i {
                            selected = true;
                        }
                    }

                    if selected {
                        if ui
                            .button(RichText::new("🔘").color(Color32::GREEN))
                            .clicked()
                        {
                            state.ui.selected_item = Selection::None;
                        }
                    } else if ui.button("🔘").clicked() {
                        state.ui.selected_item = Selection::Feature(i);
                    }

                    ui.add_space(COL_SPACING); // Less likely to accidentally delete.

                    if ui
                        .button(RichText::new("Delete 🗑").color(Color32::RED))
                        .clicked()
                    {
                        removed = Some(i);
                    }
                });

                for (key, value) in &mut feature.notes {
                    ui.horizontal(|ui| {
                        ui.label("Note:");
                        ui.add(TextEdit::singleline(key).desired_width(140.));

                        ui.label("Value:");
                        ui.add(TextEdit::singleline(value).desired_width(600.));
                    });
                }
            });

        ui.add_space(ROW_SPACING);
    }
    if let Some(rem_i) = removed {
        state.generic[state.active].features.remove(rem_i);
    }
}

pub fn feature_add_disp(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        if ui.button("➕ Add feature").clicked() {
            if state.ui.feature_add.start_posit == 0 {
                state.ui.feature_add.start_posit = 1;
            }
            if state.ui.feature_add.end_posit == 0 {
                state.ui.feature_add.end_posit = 1;
            }

            if state.ui.feature_add.start_posit > state.ui.feature_add.end_posit {
                std::mem::swap(
                    &mut state.ui.feature_add.start_posit,
                    &mut state.ui.feature_add.end_posit,
                );
            }

            state.generic[state.active].features.push(Feature {
                range: RangeIncl::new(
                    state.ui.feature_add.start_posit,
                    state.ui.feature_add.end_posit,
                ),
                feature_type: FeatureType::Generic,
                direction: FeatureDirection::None,
                label: state.ui.feature_add.label.clone(),
                color_override: None,
                notes: Default::default(),
            });
        }
    });
}

pub fn features_page(state: &mut State, ui: &mut Ui) {
    ScrollArea::vertical().show(ui, |ui| {
        feature_table(state, ui);
    });
}
