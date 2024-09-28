use core::fmt;

use eframe::egui::{Color32, ComboBox, Frame, Grid, RichText, Stroke, TextEdit, Ui, Vec2};
use strum::IntoEnumIterator;

use crate::{
    backbones::{Backbone, BackboneFilters, CloningTechnique},
    cloning::{
        find_re_candidates, make_product_tab, setup_insert_seqs, AutocloneStatus, BackboneSelected,
        CloningInsertData, Status, RE_INSERT_BUFFER,
    },
    file_io::save::load_import,
    gui::{
        find_features, lin_maps::seq_lin_disp, navigation::get_tabs, select_color_text,
        COL_SPACING, ROW_SPACING,
    },
    sequence::{seq_from_str, seq_to_str},
    util::{merge_feature_sets, RangeIncl},
    State,
};

const PASS_COLOR: Color32 = Color32::LIGHT_GREEN;
const FAIL_COLOR: Color32 = Color32::LIGHT_RED;
const NA_COLOR: Color32 = Color32::GOLD;

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

/// Draw a selector for the insert, based on loading from a file.
/// This buffer is in nucleotides, and is on either side of the insert. A buffer of 4-6 nts is ideal
/// for restriction-enzyme cloning, while no buffer is required for PCR-based cloning.
fn insert_selector(data: &mut CloningInsertData, buffer: usize, ui: &mut Ui) {
    for (i, feature) in data.features_loaded.iter().enumerate() {
        let mut border_width = 0.;
        if let Some(j) = data.feature_selected {
            if i == j {
                border_width = 1.;
            }
        }

        Frame::none()
            .stroke(Stroke::new(border_width, Color32::LIGHT_RED))
            .inner_margin(border_width)
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    if ui.button("Select").clicked {
                        data.feature_selected = Some(i);

                        // todo: Handle wraps with this for circular plasmids instead of truncating.
                        let start = if buffer + 1 < feature.range.start {
                            feature.range.start - buffer
                        } else {
                            1
                        };

                        let end_ = feature.range.end + buffer;
                        let end = if end_ + 1 < data.seq_loaded.len() {
                            end_
                        } else {
                            data.seq_loaded.len()
                        };

                        let buffered_range = RangeIncl::new(start, end);
                        if let Some(seq_this_ft) = buffered_range.index_seq(&data.seq_loaded) {
                            seq_this_ft.clone_into(&mut data.seq_insert);
                        }
                    }

                    if !feature.label.is_empty() {
                        ui.label(&feature.label);
                        ui.add_space(COL_SPACING);
                    }

                    let (r, g, b) = feature.feature_type.color();
                    ui.label(
                        RichText::new(feature.feature_type.to_string())
                            .color(Color32::from_rgb(r, g, b)),
                    );
                    ui.add_space(COL_SPACING);

                    ui.label(feature.location_descrip(data.seq_loaded.len()));
                    ui.add_space(COL_SPACING);

                    // +1 because it's inclusive.
                    ui.label(feature.location_descrip(data.seq_loaded.len()));
                });
            });
    }
}

fn insert_file_section(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        ui.label("Choose insert from:");

        let plasmid_names: &Vec<_> = &state
            .generic
            .iter()
            .map(|v| v.metadata.plasmid_name.as_str())
            .collect();

        // Add buttons for each opened tab
        for (name, i) in get_tabs(&state.path_loaded, plasmid_names, true) {
            if ui
                .button(name)
                .on_hover_text("Select an insert from this sequence")
                .clicked()
            {
                let gen = &state.generic[i];
                // This setup, including the break and variables, prevents borrow errors.
                let g = gen.features.clone();
                let s = gen.seq.clone();
                setup_insert_seqs(state, g, s);
                break;
            }
        }

        ui.add_space(COL_SPACING);
        if ui
            .button("Pick insert from file")
            .on_hover_text(
                "Choose a GenBank, PlasCAD, SnapGene, or FASTA file to \
        select an insert from. FASTA files require manual index selection.",
            )
            .clicked()
        {
            state.ui.file_dialogs.cloning_load.select_file();
        }

        ui.add_space(COL_SPACING);

        state.ui.file_dialogs.cloning_load.update(ui.ctx());

        if let Some(path) = state.ui.file_dialogs.cloning_load.take_selected() {
            if let Some(state_loaded) = load_import(&path) {
                // todo: Is there a way to do this without cloning?
                setup_insert_seqs(
                    state,
                    state_loaded.generic.features.clone(),
                    state_loaded.generic.seq.clone(),
                );
            }
        }
    });
}

fn text_from_status(status: Status) -> RichText {
    match status {
        Status::Pass => RichText::new("Pass").color(PASS_COLOR),
        Status::Fail => RichText::new("Fail").color(FAIL_COLOR),
        Status::NotApplicable => RichText::new("N/A").color(NA_COLOR),
    }
}

fn checklist(status: &AutocloneStatus, rbs_dist: Option<isize>, ui: &mut Ui) {
    ui.heading("Product checklist:");

    ui.horizontal(|ui| {
        // ui.label("Reading frame:").on_hover_text("The coding region is in-frame with respect to (todo:  RBS? Promoter?)");
        // ui.label(RichText::new("Fail").color(FAIL_COLOR));

        // todo: Handle this and missing HIS tag the same way; currently this is selectively hidden, and His is not.
        if let Some(rd) = rbs_dist {
            ui.label("Distance from RBS:").on_hover_text("Insert point is a suitabel distance (eg 5-10 nucleotides) downstream of the Ribosome Bind Site.");
            ui.label(text_from_status(status.rbs_dist));
            ui.label(format!("({rd}nt)"));
            ui.add_space(COL_SPACING);
        }
        ui.label("Downstream of promoter:").on_hover_text("Is downstream of the appropriate expression promoter.");
        ui.label(text_from_status(status.downstream_of_promoter));
        ui.add_space(COL_SPACING);

        ui.label("Upstream of terminator:").on_hover_text("Is upstream of the appropriate expression terminator.");
        ui.label(text_from_status(status.upstream_of_terminator));
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

/// A UI element that allows the user to choose which backbone to clone into.
fn backbone_selector(
    backbone_selected: &mut BackboneSelected,
    backbones: &[&Backbone],
    plasmid_name: &str,
    ui: &mut Ui,
) {
    let selected = *backbone_selected == BackboneSelected::Opened;
    if ui
        .button(select_color_text(
            &format!("This plasmid ({})", plasmid_name),
            selected,
        ))
        .clicked()
    {
        // This allows toggles.
        *backbone_selected = match backbone_selected {
            BackboneSelected::Opened => BackboneSelected::None,
            _ => BackboneSelected::Opened,
        }
    }
    ui.add_space(ROW_SPACING);

    Grid::new("0").spacing(Vec2::new(60., 6.)).show(ui, |ui| {
        for (i, backbone) in backbones.iter().enumerate() {
            let selected = match backbone_selected {
                BackboneSelected::Library(b) => *b == i,
                _ => false,
            };

            if ui
                .button(select_color_text(&backbone.name, selected))
                .clicked()
            {
                // This allows toggles.
                *backbone_selected = match backbone_selected {
                    BackboneSelected::Library(i) => BackboneSelected::None,
                    _ => BackboneSelected::Library(i),
                };
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
}

pub fn cloning_page(state: &mut State, ui: &mut Ui) {
    ui.heading("Cloning");
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

    let backbones_filtered = state.ui.backbone_filters.apply(&state.backbone_lib);

    let plasmid_name = &state.generic[state.active].metadata.plasmid_name;
    backbone_selector(
        &mut state.backbone_selected,
        &backbones_filtered,
        plasmid_name,
        ui,
    );

    let binding = Backbone::from_opened(&state.generic[state.active]);

    let bb = match state.backbone_selected {
        BackboneSelected::Library(i) => {
            if i >= state.backbone_lib.len() {
                eprintln!("Invalid index in backbone lib");
                return;
            }
            Some(&state.backbone_lib[i])
        }
        BackboneSelected::Opened => Some(&binding),
        BackboneSelected::None => None,
    };

    if let Some(backbone) = bb {
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

        let rbs_dist = backbone
            .rbs
            .map(|r| state.cloning_insert_loc as isize - r.end as isize);

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

    ui.add_space(ROW_SPACING);

    ui.horizontal(|ui| {
        ui.label("Insert location: ");
        let mut entry = state.cloning_insert_loc.to_string();
        if ui
            .add(TextEdit::singleline(&mut entry).desired_width(40.))
            .changed()
        {
            state.cloning_insert_loc = entry.parse().unwrap_or(0);
        }

        ui.add_space(COL_SPACING);
    });

    let resp_insert_editor = ui.add(
        TextEdit::multiline(&mut state.ui.cloning_insert.seq_input)
            .desired_width(ui.available_width()),
    );
    if resp_insert_editor.changed() {
        // Forces only valid NTs to be included in the string.
        state.ui.cloning_insert.seq_insert = seq_from_str(&state.ui.cloning_insert.seq_input);
        state.ui.cloning_insert.seq_input = seq_to_str(&state.ui.cloning_insert.seq_insert);
    }
}
