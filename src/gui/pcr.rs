use std::str::FromStr;

use eframe::egui::{Color32, ComboBox, Grid, RichText, TextEdit, Ui, Vec2};

use crate::{
    gui::{COL_SPACING, ROW_SPACING},
    pcr::{PolymeraseType, TempTime},
    primer::TM_TARGET,
    PcrUi, State,
};

fn temp_time_disp(tt: &TempTime, label: &str, ui: &mut Ui) {
    ui.label(&format!("{label}:"));
    ui.label(RichText::new(format!("{}°C", tt.temp)).color(Color32::LIGHT_BLUE));
    ui.label(RichText::new(format!("{}s", tt.time)).color(Color32::LIGHT_BLUE));

    ui.end_row();
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
    ui.horizontal(|ui| {
        ui.heading("PCR parameters");
        if !state.generic.primers.is_empty() {
            ui.add_space(COL_SPACING);

            if ui.button("Load from primer: ").clicked() {
                let primer = &state.generic.primers[state.ui.pcr.primer_selected]; // todo: Overflow check?

                if let Some(metrics) = &primer.volatile.metrics {
                    state.ui.pcr = PcrUi {
                        primer_tm: metrics.melting_temp,
                        product_len: primer.sequence.len(),
                        primer_selected: state.ui.pcr.primer_selected,
                        ..Default::default()
                    };
                }
                state.sync_pcr();
            }

            // Reset primer selected if an invalid one is set.
            if state.ui.pcr.primer_selected > state.generic.primers.len() {
                state.ui.pcr.primer_selected = 0;
            }

            let primer = &state.generic.primers[state.ui.pcr.primer_selected];

            ComboBox::from_id_source(0)
                .width(80.)
                .selected_text(&primer.name)
                .show_ui(ui, |ui| {
                    for (i, primer) in state.generic.primers.iter().enumerate() {
                        ui.selectable_value(&mut state.ui.pcr.primer_selected, i, &primer.name);
                    }
                });
        }
    });

    //     pub primer_tm: f32,
    //     pub product_len: usize,
    //     pub polymerase_type: PolymeraseType,
    //     pub num_cycles: u16,

    ui.horizontal(|ui| {
        // todo: Allow TM decimals?
        // Not using our helper here due to int coercing.
        ui.label("Primer TM");
        let mut entry = format!("{:.0}", state.ui.pcr.primer_tm);
        let response = ui.add(TextEdit::singleline(&mut entry).desired_width(20.));
        if response.changed() {
            state.ui.pcr.primer_tm = entry.parse().unwrap_or(TM_TARGET);
            state.sync_pcr();
        }

        // if numerical_field("Primer TM", &mut state.ui.pcr.primer_tm, TM_TARGET, ui) ||
        if numerical_field(
            "Product size (bp)",
            &mut state.ui.pcr.product_len,
            1_000,
            ui,
        ) || numerical_field("# cycles", &mut state.ui.pcr.num_cycles, 30, ui)
        {
            state.sync_pcr();
        }

        ui.label("Polymerase:");
        let prev_poly = state.ui.pcr.polymerase_type;
        ComboBox::from_id_source(1)
            .width(80.)
            .selected_text(state.ui.pcr.polymerase_type.to_str())
            .show_ui(ui, |ui| {
                ui.selectable_value(
                    &mut state.ui.pcr.polymerase_type,
                    PolymeraseType::NormalFidelity,
                    PolymeraseType::NormalFidelity.to_str(),
                );
                ui.selectable_value(
                    &mut state.ui.pcr.polymerase_type,
                    PolymeraseType::HighFidelity,
                    PolymeraseType::HighFidelity.to_str(),
                );
            });

        if state.ui.pcr.polymerase_type != prev_poly {
            state.sync_pcr();
        }
    });

    ui.add_space(ROW_SPACING);

    Grid::new(0).spacing(Vec2::new(60., 0.)).show(ui, |ui| {
        temp_time_disp(&state.pcr.initial_denaturation, "Initial denaturation", ui);
        temp_time_disp(&state.pcr.denaturation, "Denaturation", ui);
        temp_time_disp(&state.pcr.annealing, "Annealing", ui);
        temp_time_disp(&state.pcr.extension, "Extension", ui);
        temp_time_disp(&state.pcr.final_extension, "Final extension", ui);

        ui.label(format!("Number of cycles:"));
        ui.label(format!("{}", state.pcr.num_cycles));

        ui.end_row();
    });
}
