//! GUI code for the features editor and related.

use eframe::egui::{Color32, ComboBox, Frame, Painter, RichText, Sense, TextEdit, Ui, Vec2, Stroke};

use crate::{gui::{int_field, ROW_SPACING}, sequence::{
    Feature,
    FeatureDirection::{self, Forward, Reverse},
    FeatureType,
}, Color, State, Selection};

const LABEL_EDIT_WIDTH: f32 = 140.;

const COLORS: [Color; 4] = [(255, 255, 255), (255, 0, 0), (0, 255, 0), (0, 0, 255)];

/// Add a rectangle of the color for the selector.
fn color_rect(color: Color, ui: &mut Ui) {
    let color_rgb = Color32::from_rgb(color.0, color.1, color.2);

    let (rect, response) = ui.allocate_exact_size(Vec2::new(60.0, 10.0), Sense::click());
    let painter = Painter::new(ui.ctx().clone(), ui.layer_id(), rect);
    painter.rect_filled(rect, 0.0, color_rgb);
}

/// A color selector for use with feature addition and editing.
fn color_picker(val: &mut Option<Color>, id: usize, ui: &mut Ui) {
    let label_none = "Use type color";
    // let text = match val {
    //     Some((r, g, b)) => &format!("{}, {}, {}", r, g, b),
    //     None => label_none,
    // };

    color_rect(val.unwrap_or_default(), ui);

    ComboBox::from_id_source(id)
        .width(80.)
        // .selected_text(text) // todo temp
        .show_ui(ui, |ui| {
            ui.selectable_value(val, None, label_none);

            for color in COLORS {
                ui.horizontal(|ui| {
                    ui.selectable_value(val, Some(color), "");
                    color_rect(color, ui);
                });

                // ui.selectable_value(
                //     val,
                //     Some((color.0, color.1, color.2)),
                //     color.3, // todo temp. How can you show a color image?
                // );
            }
        });
}

/// A selector for use with feature addition and editing.
/// todo: Generic selector creator?
fn direction_picker(val: &mut FeatureDirection, id: usize, ui: &mut Ui) {
    ComboBox::from_id_source(id)
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
    ComboBox::from_id_source(id)
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
    for (i, feature) in state.generic.features.iter_mut().enumerate() {
        let border_color = Color32::WHITE;
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

                ui.horizontal(|ui| {

                    // todo: This may be confoudning your 0 vs 1.
                    int_field(&mut feature.index_range.0, "Start:", ui);
                    int_field(&mut feature.index_range.1, "End:", ui);

                    ui.label("Label:");
                    ui.add(TextEdit::singleline(&mut feature.label).desired_width(LABEL_EDIT_WIDTH));

                    ui.label("Type:");
                    feature_type_picker(&mut feature.feature_type, 100 + i, ui);

                    ui.label("Dir:");
                    direction_picker(&mut feature.direction, 300 + i, ui);

                    ui.label("Color:");
                    color_picker(&mut feature.color_override, 3 + i, ui);

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
                        ui.label(key); // todo!
                        // ui.add(
                        //     TextEdit::singleline(&mut note.0).desired_width(200.),
                        // );

                        ui.label("Value:");
                        ui.add(
                            TextEdit::singleline(value).desired_width(200.),
                        );
                    });
                }
            });
    }
    if let Some(rem_i) = removed {
        state.generic.features.remove(rem_i);
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

            state.generic.features.push(Feature {
                index_range: (
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
    feature_table(state, ui);
}
