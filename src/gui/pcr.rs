use std::str::FromStr;

use eframe::egui::{Color32, ComboBox, Grid, RichText, TextEdit, Ui, Vec2};

use crate::{
    gui::{COL_SPACING, ROW_SPACING},
    pcr::{PolymeraseType, TempTime},
    primer::TM_TARGET,
    PcrUi, State,
};
use crate::file_io::save::save_new_product;
use crate::gui::save::save_current_file;
use crate::primer::Primer;
use crate::util::RangeIncl;

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
    if ui
        .add(TextEdit::singleline(&mut entry).desired_width(30.))
        .changed()
    {
        *val = entry.parse().unwrap_or(default);
        changed = true;
    }

    changed
}

fn primer_dropdown(val: &mut usize, primers: &[Primer], id: usize, ui: &mut Ui) {
    // Reset primer selected if an invalid one is set.
    if *val > primers.len() {
        *val = 0;
    }

    let primer = &primers[*val];

    ComboBox::from_id_source(id)
        .width(80.)
        .selected_text(&primer.name)
        .show_ui(ui, |ui| {
            for (i, primer) in primers.iter().enumerate() {
                ui.selectable_value(val, i, &primer.name);
            }
        });
}

fn pcr_sim(state: &mut State, ui: &mut Ui) {
    let num_primers = state.generic.primers.len();

    if state.ui.pcr.primer_fwd >= num_primers || state.ui.pcr.primer_rev >= num_primers {
        state.ui.pcr.primer_fwd = 0;
        state.ui.pcr.primer_rev = 0;
    }

    ui.heading("PCR product generation");

    if num_primers >= 2 {
        ui.horizontal(|ui| {
            ui.label("Fwd:");
            primer_dropdown(&mut state.ui.pcr.primer_fwd, &state.generic.primers, 2, ui);

            ui.add_space(COL_SPACING);

            ui.label("Rev:");
            primer_dropdown(&mut state.ui.pcr.primer_rev, &state.generic.primers, 3, ui);

            ui.add_space(COL_SPACING);

            let fwd_primer = &state.generic.primers[state.ui.pcr.primer_fwd];
            let rev_primer = &state.generic.primers[state.ui.pcr.primer_rev];

            if state.ui.pcr.primer_fwd == state.ui.pcr.primer_rev {
                ui.label("Select two different primers");
            } else if  fwd_primer.volatile.matches_seq.len() == 1 && rev_primer.volatile.matches_seq.len() == 1 {
                if ui.button(RichText::new("Simulate PCR").color(Color32::GOLD)).clicked() {
                    // Save this vector; this file or quicksave instance will be turned into the PCR
                    // product.
                    save_current_file(state);

                    let fwd_primer = &state.generic.primers[state.ui.pcr.primer_fwd];
                    let rev_primer = &state.generic.primers[state.ui.pcr.primer_rev];

                    // todo: Make the the directions are correct.
                    let range_fwd = fwd_primer.volatile.matches_seq[0].1;
                    let range_rev = rev_primer.volatile.matches_seq[0].1;

                    let range_combined = RangeIncl::new(range_fwd.start, state.generic.seq.len() - range_rev.start);

                    println!("RANGE FWD: {:?}", range_fwd);
                    println!("RANGE REF: {:?}", range_rev);
                    println!("RANGE COMBINED: {:?}", range_combined);

                    if let Some(seq) = range_combined.index_seq(&state.generic.seq) {
                        state.generic.seq = seq.to_vec();
                        state.generic.features = Vec::new(); // todo: Temp
                        // state.generic.primers = Vec::new(); // todo temp.

                        state.sync_seq_related(None);

                        // todo: Remove all features not included in the range.

                        save_new_product("PCR product", state, ui);
                    }
                }
            } else {
                ui.label("There must be exactly one match for each primer");
            }
        });
    } else {
        ui.label("(Add at least 2 primers to generate a PCR product)");
    }
}

pub fn pcr_page(state: &mut State, ui: &mut Ui) {
    ui.add_space(ROW_SPACING);


    pcr_sim(state, ui);
    ui.add_space(ROW_SPACING);

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

            primer_dropdown(&mut state.ui.pcr.primer_selected, &state.generic.primers, 1, ui);
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
        ComboBox::from_id_source(10)
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

        ui.label("Number of cycles:".to_string());
        ui.label(format!("{}", state.pcr.num_cycles));

        ui.end_row();
    });
}
