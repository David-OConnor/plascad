//! GUI code related to ligation operations.

use std::collections::HashMap;

use eframe::egui::{Color32, RichText, Ui};

use crate::{
    gui::{COL_SPACING, ROW_SPACING},
    ligation::digest,
    sequence::seq_to_str,
    State,
};

pub fn ligation_page(state: &mut State, ui: &mut Ui) {
    // todo: Cache this A/R
    ui.horizontal(|ui| {
        ui.heading("Digestion and ligration");
        ui.add_space(COL_SPACING * 2.);

        ui.label("Unique cutters only:");
        ui.checkbox(&mut state.ui.re.unique_cutters_only, "");
    });

    ui.add_space(ROW_SPACING);

    // todo: Only show unique cutters A/R. And/or in the list of sites.
    // let mut res_matched = HashMap::new(); // Re, match count
    let mut res_matched = Vec::new();
    for re_match in &state.volatile[state.active].restriction_enzyme_matches {
        if re_match.lib_index + 1 >= state.restriction_enzyme_lib.len() {
            continue;
        }

        let re = &state.restriction_enzyme_lib[re_match.lib_index];

        if !res_matched.contains(&re) {
            if state.ui.re.unique_cutters_only && re_match.match_count > 1 {
                continue;
            }
            res_matched.push(re);
        }
    }

    ui.horizontal(|ui| {
        ui.heading("Restriction enzymes matched");
        ui.add_space(COL_SPACING);

        ui.label("Click to select for digestion.");
    });

    for re in res_matched {
        let selected = state.ui.re.selected.contains(&re.name);

        let color = if selected {
            Color32::LIGHT_GREEN
        } else {
            Color32::WHITE
        };

        ui.horizontal(|ui| {
            if ui.button(RichText::new(&re.name).color(color)).clicked {
                if selected {
                    for (i, name) in state.ui.re.selected.iter().enumerate() {
                        if name == &re.name {
                            state.ui.re.selected.remove(i);
                            break;
                        }
                    }
                } else {
                    state.ui.re.selected.push(re.name.clone());
                }
            }

            ui.add_space(COL_SPACING);
            ui.label(format!("{}", re.cut_depiction()));
        });

        if !state.ui.re.unique_cutters_only { // No point in displaying count for unique cutters; always 1.
             // todo: Show count.
             // ui.label(format!(": {}"));
        }
    }

    if !state.ui.re.selected.is_empty() {
        ui.add_space(ROW_SPACING);
        if ui
            .button(RichText::new("Digest").color(Color32::GOLD))
            .clicked()
        {
            state.volatile[state.active].re_digestion_products =
                digest(&state.ui.re.selected, state.get_seq());
        }
    }
}
