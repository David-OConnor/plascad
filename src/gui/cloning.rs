use core::fmt;
use std::borrow::Cow;

use eframe::egui::{
    Color32, ComboBox, Frame, Grid, RichText, ScrollArea, Stroke, TextEdit, Ui, Vec2,
};
use na_seq::{insert_into_seq, seq_from_str, seq_to_str_lower, Nucleotide};
use strum::IntoEnumIterator;

use crate::{
    backbones::{Backbone, BackboneFilters, CloningTechnique},
    cloning::{
        make_product_tab, setup_insert_seqs, BackboneSelected, CloneStatus, CloningInsertData,
        Status, RE_INSERT_BUFFER,
    },
    file_io::{save::load_import, GenericData},
    gui::{
        find_features,
        lin_maps::seq_lin_disp,
        navigation::get_tab_names,
        select_color_text,
        theme::{COLOR_ACTION, COLOR_INFO},
        COL_SPACING, ROW_SPACING,
    },
    misc_types::{Feature, FeatureType},
    state::State,
    util::{merge_feature_sets, RangeIncl},
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

    ComboBox::from_id_salt(id)
        .width(80.)
        .selected_text(text)
        .show_ui(ui, |ui| {
            ui.selectable_value(val, None, "Any");
            for variant in T::iter() {
                ui.selectable_value(val, Some(variant), variant.to_string());
            }
        });
    ui.add_space(COL_SPACING);
}

/// Descriptive text about a feature; ie for insert selection. Used to display the one selected, and in the selector.
fn draw_insert_descrip(feature: &Feature, seq_loaded: &[Nucleotide], ui: &mut Ui) {
    if !feature.label.is_empty() {
        ui.label(&feature.label);
        ui.add_space(COL_SPACING);
    }

    let (r, g, b) = feature.feature_type.color();
    ui.label(RichText::new(feature.feature_type.to_string()).color(Color32::from_rgb(r, g, b)));
    ui.add_space(COL_SPACING);

    ui.label(feature.location_descrip(seq_loaded.len()));
    ui.add_space(COL_SPACING);

    // +1 because it's inclusive.
    ui.label(feature.location_descrip(seq_loaded.len()));
}

/// Draw a selector for the insert, based on loading from a file.
/// This buffer is in nucleotides, and is on either side of the insert. A buffer of 4-6 nts is ideal
/// for restriction-enzyme cloning, while no buffer is required for PCR-based cloning.
fn insert_selector(data: &mut CloningInsertData, buffer: usize, ui: &mut Ui) -> bool {
    let mut clicked = false;

    for (i, feature) in data.features_loaded.iter().enumerate() {
        // match feature.feature_type {
        //     FeatureType::CodingRegion | FeatureType::Generic | FeatureType::Gene => (),
        //     _ => continue,
        // }

        let mut selected = false;
        if let Some(j) = data.feature_selected {
            if i == j {
                selected = true;
            }
        }

        let (border_width, btn_color) = if selected {
            (1., Color32::GREEN)
        } else {
            (0., Color32::WHITE)
        };

        Frame::none()
            .stroke(Stroke::new(border_width, Color32::LIGHT_RED))
            .inner_margin(border_width)
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    if ui
                        .button(RichText::new("Select").color(btn_color))
                        .clicked()
                    {
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
                        clicked = true;
                    }

                    draw_insert_descrip(feature, &data.seq_loaded, ui);
                });
            });
    }
    clicked
}

/// Choose from tabs to select an insert from.
fn insert_tab_selection(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        ui.label("Choose insert from:");

        let plasmid_names: &Vec<_> = &state
            .generic
            .iter()
            .map(|v| v.metadata.plasmid_name.as_str())
            .collect();

        // Add buttons for each opened tab
        for (name, i) in get_tab_names(&state.tabs_open, plasmid_names, true) {
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

                state.ui.cloning_insert.show_insert_picker = true;

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
            state.ui.file_dialogs.cloning_load.pick_file();
            state.ui.cloning_insert.show_insert_picker = true;
        }

        ui.add_space(COL_SPACING);

        state.ui.file_dialogs.cloning_load.update(ui.ctx());

        if let Some(path) = state.ui.file_dialogs.cloning_load.take_picked() {
            if let Some(state_loaded) = load_import(&path) {
                // todo: Is there a way to do this without cloning?
                setup_insert_seqs(
                    state,
                    state_loaded.generic.features.clone(),
                    state_loaded.generic.seq.clone(),
                );
            }
        }

        if !state.ui.cloning_insert.features_loaded.is_empty() {
            let hide_text = if state.ui.cloning_insert.show_insert_picker {
                "Hide inserts"
            } else {
                "Show inserts"
            };
            ui.add_space(COL_SPACING);

            if ui.button(hide_text).clicked() {
                state.ui.cloning_insert.show_insert_picker =
                    !state.ui.cloning_insert.show_insert_picker
            }
        }
    });

    // ui.add_space(ROW_SPACING);

    ui.horizontal(|ui| {
        // A short summary of the selected feature; useful if the picker is hidden.
        for (i, feature) in state.ui.cloning_insert.features_loaded.iter().enumerate() {
            let mut border_width = 0.;
            if let Some(j) = state.ui.cloning_insert.feature_selected {
                if i == j {
                    ui.label(RichText::new("Insert selected:").color(COLOR_INFO));
                    draw_insert_descrip(feature, &state.ui.cloning_insert.seq_loaded, ui);
                }
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

fn checklist(status: &CloneStatus, rbs_dist: Option<isize>, ui: &mut Ui) {
    ui.heading("Product checklist:");

    ui.horizontal(|ui| {
        // ui.label("Reading frame:").on_hover_text("The coding region is in-frame with respect to (todo:  RBS? Promoter?)");
        // ui.label(RichText::new("Fail").color(FAIL_COLOR));

        ui.label("Distance from RBS:").on_hover_text("Insert point is a suitabel distance (eg 5-10 nucleotides) downstream of the Ribosome Bind Site.");
        ui.label(text_from_status(status.rbs_dist));
        if let Some(rd) = rbs_dist {
            ui.label(format!("({rd}nt)"));
        }
        ui.add_space(COL_SPACING);

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
    data: &GenericData,
    bb_cache: &mut Option<Backbone>,
    ui: &mut Ui,
) -> bool {
    let mut changed = false;

    ui.horizontal(|ui| {
        if ui
            .button(select_color_text(
                &format!("This plasmid ({})", plasmid_name),
                *backbone_selected == BackboneSelected::Opened,
            ))
            .clicked()
        {
            // Cache this, since we don't have it in the library to reference.
            *bb_cache = Some(Backbone::from_opened(data));
            *backbone_selected = BackboneSelected::Opened;
            changed = true;
        }
        ui.add_space(COL_SPACING);

        if backbones.is_empty() {
            eprintln!("Error: Empty backbone library");
            return;
        }

        if ui
            .button(select_color_text(
                "Library:",
                *backbone_selected != BackboneSelected::Opened,
            ))
            .clicked()
        {
            *backbone_selected = BackboneSelected::Library(0);
            changed = true
        }
        ui.add_space(COL_SPACING);

        let bb_selected = match backbone_selected {
            BackboneSelected::Opened => backbones[0],
            BackboneSelected::Library(i) => backbones[*i],
        };

        let bb_prev = &backbone_selected.clone(); // todo: Don't like this clone.
        ComboBox::from_id_salt(1000)
            .width(80.)
            .selected_text(&bb_selected.name)
            .show_ui(ui, |ui| {
                for (i, backbone) in backbones.iter().enumerate() {
                    ui.selectable_value(
                        backbone_selected,
                        BackboneSelected::Library(i),
                        &backbone.name,
                    );
                }
            });

        if *backbone_selected != *bb_prev {
            changed = true;
        }

        ui.add_space(COL_SPACING);

        if let Some(addgene_url) = bb_selected.addgene_url() {
            if ui.button("View on AddGene").clicked() {
                if let Err(e) = webbrowser::open(&addgene_url) {
                    eprintln!("Failed to open the web browser: {:?}", e);
                }
            }
        }
    });

    ui.add_space(ROW_SPACING);

    changed
}

pub fn cloning_page(state: &mut State, ui: &mut Ui) {
    ScrollArea::vertical().id_salt(100).show(ui, |ui| {
        let mut sync = false;

        ui.heading("Cloning (Currently supports PCR-based cloning only");
        //     ui.label("For a given insert, automatically select a backbone, and either restriction enzymes, or PCR primers to use\
        // to clone the insert into the backbone.");

        ui.add_space(ROW_SPACING);

        ui.horizontal(|ui| {
            ui.label(
                "Remove stop coding prior to tags. (Useful if there is a tag on the backbone)",
            )
            .on_hover_text(
                "If there is a stop coding at the end of the insert, and a His or similar tag, \
                remove the codon, so the tag is coded for",
            );
            if ui
                .checkbox(&mut state.cloning.remove_stop_codons, "")
                .changed()
            {
                sync = true;
            }
        });
        ui.add_space(ROW_SPACING);

        // todo: DRY with below getting the backbone, but we have a borrow error when moving that up.
        let data_vec = match state.cloning.backbone_selected {
            BackboneSelected::Library(i) => {
                if i >= state.backbone_lib.len() {
                    eprintln!("Invalid index in backbone lib");
                    None
                } else {
                    Some(&state.backbone_lib[i].data)
                }
            }
            BackboneSelected::Opened => Some(&state.generic[state.active]),
        };

        // Draw the linear map regardless of if there's a vector (Empty map otherwise). This prevents
        // a layout shift when selecting.
        let data_vec_ = if let Some(data_) = data_vec {
            data_
        } else {
            &Default::default()
        };

        // A minimap for the vector.
        seq_lin_disp(
            data_vec_,
            true,
            state.ui.selected_item,
            &state.ui.re.res_selected,
            Some(state.cloning.insert_loc),
            &state.ui,
            &state.cloning.re_matches_vec_common,
            &state.restriction_enzyme_lib,
            ui,
        );

        ui.add_space(ROW_SPACING);

        // A minimap for the insert
        if let Some(data) = &state.cloning.data_insert {
            seq_lin_disp(
                data,
                true,
                state.ui.selected_item,
                &state.ui.re.res_selected,
                None,
                &state.ui,
                &state.cloning.re_matches_insert_common,
                &state.restriction_enzyme_lib,
                ui,
            );
        }

        ui.add_space(ROW_SPACING);

        insert_tab_selection(state, ui);
        ui.add_space(ROW_SPACING);

        if state.ui.cloning_insert.show_insert_picker {
            let insert_just_picked =
                insert_selector(&mut state.ui.cloning_insert, RE_INSERT_BUFFER, ui);

            if insert_just_picked {
                sync = true;
            }
            ui.add_space(ROW_SPACING);
        }

        ui.heading("Backbones (vectors)");
        backbone_filters(&mut state.ui.backbone_filters, ui);
        ui.add_space(ROW_SPACING);

        // todo: Cache this?
        let backbones_filtered = state.ui.backbone_filters.apply(&state.backbone_lib);

        let plasmid_name = &state.generic[state.active].metadata.plasmid_name;
        let backbone_just_picked = backbone_selector(
            &mut state.cloning.backbone_selected,
            &backbones_filtered,
            plasmid_name,
            &state.generic[state.active],
            &mut state.cloning.backbone,
            ui,
        );

        if backbone_just_picked {
            sync = true;
        }

        // todo: This is DRY with get_backbone due to borrow error.
        // let backbone = state.cloning.get_backbone(&state.backbone_lib);
        let backbone = match state.cloning.backbone_selected {
            BackboneSelected::Library(i) => {
                if i >= state.backbone_lib.len() {
                    eprintln!("Invalid index in backbone lib");
                    None
                } else {
                    Some(&state.backbone_lib[i])
                }
            }
            BackboneSelected::Opened => match state.cloning.backbone.as_ref() {
                Some(bb) => Some(bb),
                None => None,
            },
        };

        // These variables prevent borrow errors on backbone.
        let mut clone_initiated = false;

        if let Some(backbone) = backbone {
            let rbs_dist = backbone
                .rbs
                .map(|r| state.cloning.insert_loc as isize - r.end as isize);

            ui.add_space(ROW_SPACING);
            ui.label("Restriction enzymes:");
            if state.cloning.res_common.is_empty() {
                ui.label("(None)");
            }

            for candidate in &state.cloning.res_common {
                ui.label(&candidate.name);
            }

            ui.add_space(ROW_SPACING);

            if ui
                .button(
                    RichText::new("Auto set insert location (PCR; expression)").color(COLOR_ACTION),
                )
                .clicked()
            {
                if let Some(insert_loc) = backbone.insert_loc(CloningTechnique::Pcr) {
                    state.cloning.insert_loc = insert_loc;
                }
                sync = true;
            }

            ui.horizontal(|ui| {
                ui.label("Insert location:");
                ui.label(RichText::new(format!("{}", state.cloning.insert_loc)).color(COLOR_INFO));
            });
            ui.add_space(ROW_SPACING);

            // todo: Only if there is a result
            if true {
                ui.add_space(ROW_SPACING);
                checklist(&state.cloning.status, rbs_dist, ui);

                ui.add_space(ROW_SPACING);

                if ui
                    .button(RichText::new("Clone (PCR)").color(COLOR_ACTION))
                    .clicked()
                {
                    clone_initiated = true;
                }
            }

            ui.add_space(ROW_SPACING);

            ui.horizontal(|ui| {
                ui.label("Insert location: ");
                let mut entry = state.cloning.insert_loc.to_string();
                if ui
                    .add(TextEdit::singleline(&mut entry).desired_width(40.))
                    .changed()
                {
                    state.cloning.insert_loc = entry.parse().unwrap_or(0);
                    sync = true;
                }

                ui.add_space(COL_SPACING);
            });

            if clone_initiated {
                make_product_tab(state, Some(backbone.data.clone()));
                // Annotate the vector, for now at least.
                let features = find_features(&state.get_seq());
                // We assume the product has been made active.
                merge_feature_sets(&mut state.generic[state.active].features, &features)
            }
        }

        let resp_insert_editor = ui.add(
            TextEdit::multiline(&mut state.ui.cloning_insert.seq_input)
                .desired_width(ui.available_width()),
        );
        if resp_insert_editor.changed() {
            // Forces only valid NTs to be included in the string.
            state.ui.cloning_insert.seq_insert = seq_from_str(&state.ui.cloning_insert.seq_input);
            state.ui.cloning_insert.seq_input =
                seq_to_str_lower(&state.ui.cloning_insert.seq_insert);

            sync = true;
        }

        if sync {
            state.cloning.sync(
                &mut state.ui.cloning_insert.seq_insert,
                &state.backbone_lib,
                &state.restriction_enzyme_lib,
            );
        }
    });
}
