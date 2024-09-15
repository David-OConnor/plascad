use core::fmt;

use eframe::egui::{ComboBox, Ui};
use strum::IntoEnumIterator;
use strum_macros::EnumIter;

use crate::{
    autocloning::{find_re_candidates, RE_INSERT_BUFFER},
    backbones,
    backbones::{BackboneFilters, ExpressionHost},
    gui::{
        cloning::{insert_file_section, insert_selector},
        lin_maps::seq_lin_disp,
        select_color_text, COL_SPACING, ROW_SPACING,
    },
    State,
};

fn filter_selector<T: fmt::Display + PartialEq + Copy + IntoEnumIterator>(
    name: &str,
    val: &mut Option<T>,
    id: u32,
    ui: &mut Ui,
) {
    ui.label(name);

    let text = match val {
        Some(v) => v.to_string(),
        None => "Any".to_string(),
    };

    ComboBox::from_id_source(id)
        .width(80.)
        .selected_text(text)
        .show_ui(ui, |ui| {
            ui.selectable_value(val, None, "Any");
            // todo: Impl Iter on the enums
            for variant in T::iter() {
                ui.selectable_value(val, Some(variant), variant.to_string());
            }
        });
    ui.add_space(COL_SPACING);
}

fn backbone_filters(filters: &mut BackboneFilters, ui: &mut Ui) {
    // todo: Allow selecting multiple options.

    // todo: HElper fn for these
    ui.horizontal(|ui| {
        filter_selector("Host", &mut filters.host, 0, ui);
        filter_selector("Ab resistance", &mut filters.antibiotic_resistance, 1, ui);
        filter_selector("Expression", &mut filters.expression_system, 2, ui);
        filter_selector("Copy", &mut filters.copy_number, 3, ui);

        ui.label("His tagged:");
        ui.checkbox(&mut filters.his_tagged, "");
        ui.add_space(COL_SPACING);
    });
}

pub fn autocloning_page(state: &mut State, ui: &mut Ui) {
    ui.heading("Assisted cloning");
    ui.label("For a given insert, automatically select a backbone, and either restriction enzymes, or PCR primers to use\
    to clone the insert into the backbone.");

    ui.add_space(ROW_SPACING);

    // todo: Replace this with something more appropriate; ie a sub-call of it etc that shows teh backbone
    // todo or backbone with insert.
    seq_lin_disp(state, ui, true, state.active, &state.ui.re.res_selected);
    ui.add_space(ROW_SPACING);

    insert_file_section(state, ui);
    ui.add_space(ROW_SPACING);

    insert_selector(&mut state.ui.cloning_insert, RE_INSERT_BUFFER, ui);
    ui.add_space(ROW_SPACING);

    ui.heading("Backbones (vectors)");
    backbone_filters(&mut state.ui.backbone_filters, ui);
    ui.add_space(ROW_SPACING);

    let mut backbones_filtered = state.ui.backbone_filters.apply(&state.backbone_lib);

    for (i, backbone) in backbones_filtered.iter().enumerate() {
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

            if let Some(addgene_url) = backbone.addgene_url() {
                if ui.button("View on AddGene").clicked() {
                    if let Err(e) = webbrowser::open(&addgene_url) {
                        eprintln!("Failed to open the web browser: {:?}", e);
                    }
                }
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
            &state.ui.cloning_insert.seq_insert,
            &state.restriction_enzyme_lib,
            &state.volatile,
        );

        ui.add_space(ROW_SPACING);
        ui.label("Restriction enzymes:");
        if state.cloning_res_matched.is_empty() {
            ui.label("(None)");
        }

        for candidate in &state.cloning_res_matched {
            ui.label(&candidate.name);
        }
    }
}
