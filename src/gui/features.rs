//! GUI code for the features editor and related.

use eframe::{
    egui::{
        pos2, Align2, Color32, ComboBox, FontFamily, FontId, Pos2, RichText, Shape, Stroke,
        TextEdit, Ui,
    },
    epaint::PathShape,
};

use crate::{
    gui::int_field,
    sequence::{
        Feature,
        FeatureDirection::{self, Forward, Reverse},
        FeatureType,
    },
    Color, State,
};

const LABEL_EDIT_WIDTH: f32 = 140.;

// todo: The str label here is temp; ideally have the picker show the color itself.
const COLORS: [(u8, u8, u8, &str); 4] = [
    (255, 255, 255, "White"),
    (255, 0, 0, "Red"),
    (0, 255, 0, "Blue"),
    (0, 0, 255, "Green"),
];

/// A color selector for use with feature addition and editing.
fn color_picker(val: &mut Option<Color>, id: usize, ui: &mut Ui) {
    let label_none = "Use type color";
    let text = match val {
        Some((r, g, b)) => &format!("{}, {}, {}", r, g, b),
        None => label_none,
    };
    ComboBox::from_id_source(id)
        .width(80.)
        .selected_text(text) // todo temp
        .show_ui(ui, |ui| {
            ui.selectable_value(val, None, label_none);
            for color in COLORS {
                ui.selectable_value(
                    val,
                    Some((color.0, color.1, color.2)),
                    color.3, // todo temp. How can you show a color image?
                );
            }
        });
}

/// A selector for use with feature addition and editing.
/// todo: Generic selector creator?
fn feature_type_picker(val: &mut FeatureType, id: usize, ui: &mut Ui) {
    ComboBox::from_id_source(id)
        .width(140.)
        .selected_text(val.to_string())
        .show_ui(ui, |ui| {
            for feature_type in [
                FeatureType::Generic,
                FeatureType::Gene,
                FeatureType::Ori,
                // FeatureType::RnaPolyBindSite,
                FeatureType::RibosomeBindSite,
                FeatureType::AntibioticResistance,
                FeatureType::CodingRegion,
            ] {
                ui.selectable_value(val, feature_type, feature_type.to_string());
            }
        });
}

pub fn feature_table(features: &mut Vec<Feature>, ui: &mut Ui) {
    let mut removed = None;
    for (i, feature) in features.iter_mut().enumerate() {
        ui.horizontal(|ui| {
            // todo: This may be confoudning your 0 vs 1.
            int_field(&mut feature.index_range.0, "Start:", ui);
            int_field(&mut feature.index_range.1, "End:", ui);

            ui.label("Label:");
            ui.add(TextEdit::singleline(&mut feature.label).desired_width(LABEL_EDIT_WIDTH));

            ui.label("Type:");
            feature_type_picker(&mut feature.feature_type, 100 + i, ui);

            ui.label("Color:");
            color_picker(&mut feature.color_override, 3 + i, ui);

            if ui
                .button(RichText::new("Delete 🗑").color(Color32::RED))
                .clicked()
            {
                removed = Some(i);
            }
        });
    }
    if let Some(rem_i) = removed {
        features.remove(rem_i);
    }
}

pub fn feature_add_disp(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        ui.heading("Add feature: ");

        int_field(&mut state.ui.feature_add.start_posit, "Start:", ui);
        int_field(&mut state.ui.feature_add.end_posit, "End:", ui);

        ui.label("Label:");
        ui.add(
            TextEdit::singleline(&mut state.ui.feature_add.label).desired_width(LABEL_EDIT_WIDTH),
        );

        ui.label("Type:");
        feature_type_picker(&mut state.ui.feature_add.feature_type, 200, ui);

        ui.label("Color:");
        color_picker(&mut state.ui.feature_add.color, 2, ui);

        if ui.button("➕ Add").clicked() {
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

            state.features.push(Feature {
                index_range: (
                    state.ui.feature_add.start_posit,
                    state.ui.feature_add.end_posit,
                ),
                feature_type: FeatureType::Generic,
                direction: FeatureDirection::None,
                label: state.ui.feature_add.label.clone(),
                color_override: None,
            });
        }
    });
}
