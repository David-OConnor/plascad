use eframe::egui::{RichText, Ui};

use crate::{
    autocloning::find_re_candidates,
    gui::{
        cloning::{insert_file_section, insert_selector},
        select_color_text, ROW_SPACING,
    },
    State,
};

pub fn autocloning_page(state: &mut State, ui: &mut Ui) {
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

    if let Some(bb_i) = state.backbone_selected {
        // todo: Cache any calcs you are currently doing here.
        if bb_i >= state.backbone_lib.len() {
            eprintln!("Invalid index in backbone lib");
            return;
        }
        let backbone = &state.backbone_lib[bb_i];

        // todo: Move where this is run so it's only run when appropriate! May need some additiosn to the state sync fns as well.
        state.cloning_res_matched = find_re_candidates(
            &backbone,
            // &state.ui.cloning_insert.seq_insert,
            &state.restriction_enzyme_lib,
            &state.volatile,
        );

        ui.label("Restriction enzymes");
        for candidate in &state.cloning_res_matched {

        }
    }
}
