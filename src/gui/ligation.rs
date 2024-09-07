//! GUI code related to ligation operations.

use std::collections::HashMap;

use eframe::egui::Ui;

use crate::{gui::COL_SPACING, State};

pub fn ligation_page(state: &mut State, ui: &mut Ui) {
    // todo: Cache this A/R
    ui.horizontal(|ui| {
        ui.heading("Digestion and ligration");
        ui.add_space(COL_SPACING * 2.);

        ui.label("Unique cutters only:");
        ui.checkbox(&mut state.ui.re.unique_cutters_only, "");
    });

    // todo: Only show unique cutters A/R. And/or in the list of sites.
    let mut res_matched = HashMap::new(); // Re, match count
    for re_match in &state.volatile.restriction_enzyme_matches {
        if re_match.lib_index + 1 >= state.restriction_enzyme_lib.len() {
            continue;
        }

        let re = &state.restriction_enzyme_lib[re_match.lib_index];

        if res_matched.contains_key(&re) {
            *res_matched.get_mut(&re).unwrap() += 1;
        } else {
            res_matched.insert(re, 0);
        }
    }

    ui.heading("Restriction enzymes matched:");
    ui.horizontal(|ui| {
        for (re, count) in res_matched {
            if state.ui.re.unique_cutters_only && count > 1 {
                continue;
            }
            ui.label(&re.name);
        }
    });
}
