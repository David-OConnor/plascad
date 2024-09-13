use eframe::egui::RichText;

use crate::{
    gui::{
        cloning::{insert_file_section, insert_selector},
        select_color_text, ROW_SPACING,
    },
    State,
};

pub fn autocloning_page(state: &mut State, ui: &mut eframe::egui::Ui) {
    ui.heading("Assisted cloning");
    ui.label("For a given insert, automatically select a backbone, and either restriction enzymes, or PCR primers to use\
    to clone the insert into the backbone.");

    ui.add_space(ROW_SPACING);

    insert_file_section(state, ui);
    ui.add_space(ROW_SPACING);

    insert_selector(&mut state.ui.cloning_insert, ui);
    ui.add_space(ROW_SPACING);

    ui.heading("Backbones (vectors)");

    for (i, backbone) in state.backbone_lib.iter().enumerate() {
        ui.horizontal(|ui| {
            let selected = match state.backbone_selected {
                Some(b) => b == i,
                None => false,
            };

            // todo: Fn for this button selecting/green
            if ui
                .button(select_color_text(&backbone.name, selected))
                .clicked()
            {
                // todo: Unselect if selected
                state.backbone_selected = Some(i);
            }
        });
    }
}
