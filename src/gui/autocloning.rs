use core::fmt;

use eframe::egui::{Color32, ComboBox, Grid, RichText, Ui, Vec2};
use strum::IntoEnumIterator;

use crate::{gui::find_features, util::merge_feature_sets};

const PASS_COLOR: Color32 = Color32::LIGHT_GREEN;
const FAIL_COLOR: Color32 = Color32::LIGHT_RED;

use crate::{
    autocloning::{find_re_candidates, AutocloneStatus, Status, RE_INSERT_BUFFER},
    backbones::{BackboneFilters, CloningTechnique},
    cloning::make_product_tab,
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

fn text_from_status(status: Status) -> RichText {
    match status {
        Status::Pass => RichText::new("Pass").color(PASS_COLOR),
        Status::Fail => RichText::new("Fail").color(FAIL_COLOR),
    }
}

fn checklist(status: &AutocloneStatus, rbs_dist: isize, ui: &mut Ui) {
    ui.heading("Product checklist:");

    ui.horizontal(|ui| {
        // ui.label("Reading frame:").on_hover_text("The coding region is in-frame with respect to (todo:  RBS? Promoter?)");
        // ui.label(RichText::new("Fail").color(FAIL_COLOR));

        ui.label("Distance from RBS:").on_hover_text("Insert point is a suitabel distance (eg 5-10 nucleotides) downstream of the Ribosome Bind Site.");
        ui.label(text_from_status(status.rbs_dist));
        ui.label(format!("({rbs_dist}nt)"));
        ui.add_space(COL_SPACING);

        ui.label("Downstream of promoter:").on_hover_text("Is downstream of the appropriate expression promoter.");
        ui.label(text_from_status(status.downstream_of_promoter));
        ui.add_space(COL_SPACING);

        ui.label("Direction:");
        ui.label(text_from_status(status.direction));
        ui.add_space(COL_SPACING);

        ui.label("In frame with His tag:");
        ui.label(text_from_status(status.tag_frame));
        ui.add_space(COL_SPACING);

        // ui.label("Primer quality:");
        // ui.label(RichText::new("Fail").color(FAIL_COLOR));
        // ui.add_space(COL_SPACING);
        //
        // ui.label("RE distance:");
        // ui.label(RichText::new("Fail").color(FAIL_COLOR));
        // ui.add_space(COL_SPACING);
    });
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
    ui.heading("Assisted cloning (Work in progress)");
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
    Grid::new("0").spacing(Vec2::new(60., 6.)).show(ui, |ui| {
        for (i, backbone) in backbones_filtered.iter().enumerate() {
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
            ui.end_row();
        }
    });

    if let Some(bb_i) = state.backbone_selected {
        if bb_i >= state.backbone_lib.len() {
            eprintln!("Invalid index in backbone lib");
            return;
        }
        let backbone = &state.backbone_lib[bb_i];

        // todo: Cache all relevant calcs you are currently doing here! For example, only when you change vector,
        // todo: or when the sequence changes. (Eg with state sync fns)
        state.cloning_res_matched = find_re_candidates(
            &backbone,
            &state.ui.cloning_insert.seq_insert,
            &state.restriction_enzyme_lib,
            &state.volatile,
        );

        if let Some(insert_loc) = backbone.insert_loc(CloningTechnique::Pcr) {
            state.cloning_insert_loc = insert_loc;
        }

        state.autoclone_status = AutocloneStatus::new(
            &backbone,
            state.cloning_insert_loc,
            state.ui.cloning_insert.seq_insert.len(),
        );

        let rbs_dist = state.cloning_insert_loc as isize - backbone.rbs.end as isize;

        // todo: End calcs to cache.

        ui.add_space(ROW_SPACING);
        ui.label("Restriction enzymes:");
        if state.cloning_res_matched.is_empty() {
            ui.label("(None)");
        }

        for candidate in &state.cloning_res_matched {
            ui.label(&candidate.name);
        }

        ui.add_space(ROW_SPACING);

        if let Some(insert_loc) = backbone.insert_loc(CloningTechnique::Pcr) {
            ui.horizontal(|ui| {
                ui.label("Insert location (PCR):");
                ui.label(RichText::new(format!("{insert_loc}")).color(Color32::LIGHT_BLUE));
            });
        }
        ui.add_space(ROW_SPACING);

        // todo: Only if there is a result
        if true {
            ui.add_space(ROW_SPACING);
            checklist(&state.autoclone_status, rbs_dist, ui);

            ui.add_space(ROW_SPACING);

            if ui
                .button(RichText::new("Clone (PCR)").color(Color32::GOLD))
                .clicked()
            {
                make_product_tab(state, Some(backbone.make_generic_data()));
                // Annotate the vector, for now at least.
                let features = find_features(&state.get_seq());
                // We assume the product has been made active.
                merge_feature_sets(&mut state.generic[state.active].features, &features)
            }
        }
    }
}
