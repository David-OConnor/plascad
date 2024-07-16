use std::str::FromStr;
use eframe::egui::{ComboBox, TextEdit, Ui};
use crate::gui::{COL_SPACING, ROW_SPACING};
use crate::pcr::{PcrParams, PolymeraseType, TempTime};
use crate::primer::TM_TARGET;
use crate::State;

fn temp_time_disp(tt: &TempTime, label: &str, ui: &mut Ui) {
    ui.horizontal(|ui| {
        ui.label(&format!("{label}:"));
        ui.add_space(COL_SPACING);

        ui.label(&format!("{}°C", tt.temp));
        ui.add_space(COL_SPACING);

        ui.label(&format!("{}s", tt.time));
        ui.add_space(COL_SPACING);

    });
}

fn numerical_field<T>(label: &str, val: &mut T, default: T, ui: &mut Ui) -> bool
where
    T: ToString + FromStr + Copy,
    <T as FromStr>::Err: std::fmt::Debug,
{
    let mut changed = false;

    ui.label(&format!("{}: ", label));
    let mut entry = val.to_string();
    let response = ui.add(TextEdit::singleline(&mut entry).desired_width(30.));
    if response.changed() {
        *val = entry.parse().unwrap_or(default);
        changed = true;
    }

    changed
}

pub fn pcr_page(state: &mut State, ui: &mut Ui) {
    ui.heading("PCR parameters");

    //     pub primer_tm: f32,
    //     pub product_len: usize,
    //     pub polymerase_type: PolymeraseType,
    //     pub num_cycles: u16,

    ui.horizontal(|ui| {
        let mut changed = state.sync_pcr();
        if numerical_field("Primer TM", &mut state.ui.pcr.primer_tm, TM_TARGET, ui) ||
            numerical_field("Product size (bp)", &mut state.ui.pcr.product_len, 0, ui) ||
            numerical_field("# cycles", &mut state.ui.pcr.num_cycles, 30, ui)
        {
            state.sync_pcr();
        }


        ui.label("Polymerase:");
        let prev_poly = state.ui.pcr.polymerase_type;
        ComboBox::from_id_source(0)
            .width(80.)
            .selected_text(state.ui.pcr.polymerase_type.to_str())
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut state.ui.pcr.polymerase_type, PolymeraseType::NormalFidelity, PolymeraseType::NormalFidelity.to_str());
                ui.selectable_value(&mut state.ui.pcr.polymerase_type, PolymeraseType::HighFidelity, PolymeraseType::HighFidelity.to_str());
            });

        if state.ui.pcr.polymerase_type != prev_poly {
            state.sync_pcr();
        }

        // ui.label("Primer TM");
        // let mut val = state.ui.pcr.primer_tm;
        // let mut entry = val.to_string();
        // let response = ui.add(TextEdit::singleline(&mut entry).desired_width(40.));
        // if response.changed() {
        //     state.ui.pcr.primer_tm = entry.parse().unwrap_or(TM_TARGET);
        //     state.sync_pcr();
        // }
        // ui.add_space(COL_SPACING);
        //
        // ui.label();
        // let mut val = state.ui.pcr.product_len;
        // let mut entry = val.to_string();
        // let response = ui.add(TextEdit::singleline(&mut entry).desired_width(40.));
        // if response.changed() {
        //     state.ui.pcr.product_len = entry.parse().unwrap_or(1_000);
        //     state.sync_pcr();
        // }

    });

    ui.add_space(ROW_SPACING);



    temp_time_disp(&state.pcr.initial_denaturation, "Initial denaturation", ui);
    temp_time_disp(&state.pcr.denaturation, "Denaturation", ui);
    temp_time_disp(&state.pcr.annealing, "Annealing", ui);
    temp_time_disp(&state.pcr.extension, "Extension", ui);
    temp_time_disp(&state.pcr.final_extension, "Final extension", ui);

    ui.horizontal(|ui| {
        ui.label(&format!("Number of cycles:"));
        ui.add_space(COL_SPACING);

        ui.label(&format!("{}", state.pcr.num_cycles));
        ui.add_space(COL_SPACING);
    });
}
