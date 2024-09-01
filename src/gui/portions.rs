//! UI page for mixing portions. (growth media, stock solutions etc)

use eframe::{
    egui,
    egui::{Color32, ComboBox, RichText, TextEdit, Ui},
};

use crate::{
    gui::{COL_SPACING, ROW_SPACING},
    portions::{
        media_prep, MediaPrepInput, PlateSize, PortionsState, Reagent, ReagentPrep, ReagentType,
        Solution,
    },
};
// todo: Make a non-gui portions module once this becomes unweildy.

// todo: Store solutions. Save to file, and mix solutions from other solutions.

const DEFAULT_REAGENT_MOLARITY: f32 = 1.;
const DEFAULT_TOTAL_VOLUME: f32 = 1.; // L

fn solutions_disp(portions: &mut PortionsState, ui: &mut Ui) {
    let mut sol_removed = None;
    for (i, solution) in portions.solutions.iter_mut().enumerate() {
        ui.horizontal(|ui| {
            // ui.heading(RichText::new(&solution.name).color(Color32::LIGHT_BLUE));
            ui.add(TextEdit::singleline(&mut solution.name).desired_width(200.));

            ui.label("Total volume (mL):");

            // todo: Float field, or can we do ml?
            let mut val_ml = ((solution.total_volume * 1_000.) as u32).to_string();
            let response = ui.add(TextEdit::singleline(&mut val_ml).desired_width(40.));

            if response.changed() {
                let val_int: u32 = val_ml.parse().unwrap_or_default();
                solution.total_volume = val_int as f32 / 1_000.;
                solution.calc_amounts();
            }

            if ui
                .button(RichText::new("➕ Add reagent").color(Color32::GOLD))
                .clicked()
            {
                solution.reagents.push(Reagent::default());
            }

            if ui
                .button(RichText::new("Delete solution 🗑").color(Color32::RED))
                .clicked()
            {
                sol_removed = Some(i);
            }
        });

        ui.add_space(ROW_SPACING / 2.);

        let mut reagent_removed = None;

        for (j, reagent) in solution.reagents.iter_mut().enumerate() {
            ui.horizontal(|ui| {
                let type_prev = reagent.type_.clone();
                ComboBox::from_id_source(100 + i * 100 + j)
                    .width(140.)
                    .selected_text(reagent.type_.to_string())
                    .show_ui(ui, |ui| {
                        for type_ in [
                            ReagentType::Custom(0.),
                            ReagentType::SodiumChloride,
                            ReagentType::TrisHcl,
                            ReagentType::Iptg,
                            ReagentType::SodiumPhosphateMonobasic,
                            ReagentType::SodiumPhosphateDibasic,
                            ReagentType::SodiumPhosphateDibasicHeptahydrate,
                            ReagentType::PotassiumPhosphateMonobasic,
                            ReagentType::PotassiumPhosphateDibasic,
                            ReagentType::Imidazole,
                            ReagentType::Lysozyme,
                            ReagentType::Mes,
                            ReagentType::Bes,
                            ReagentType::Tes,
                            ReagentType::CitricAcid,
                            ReagentType::Edta,
                            ReagentType::HydrochloricAcid,
                            ReagentType::SodiumHydroxide,
                        ] {
                            ui.selectable_value(&mut reagent.type_, type_, type_.to_string());
                        }

                        // ui.selectable_value(val, Reagent::Custom(_), Reagent::Custom(_).to_string());
                    });

                if let ReagentType::Custom(weight) = &mut reagent.type_ {
                    let mut val = format!("{:.2}", weight);
                    if ui
                        .add(TextEdit::singleline(&mut val).desired_width(40.))
                        .changed
                    {
                        // let molarity_int: u32 = val.parse().unwrap_or_default();
                        *weight = val.parse().unwrap_or_default();
                        // *v = molarity_int as f32 / 1_000.;
                        reagent.calc_amount(solution.total_volume);
                    }
                    ui.label("g/mol");
                } else {
                    ui.add_sized(
                        [80.0, 20.0],
                        egui::Label::new(format!("{:.2} g/mol", reagent.type_.weight())),
                    );
                }

                let prep_prev = reagent.prep.clone();
                ComboBox::from_id_source(2000 + i * 100 + j)
                    .width(80.)
                    .selected_text(reagent.prep.to_string())
                    .show_ui(ui, |ui| {
                        for prep in [
                            ReagentPrep::Mass,
                            ReagentPrep::Volume(DEFAULT_REAGENT_MOLARITY),
                        ] {
                            ui.selectable_value(&mut reagent.prep, prep, prep.to_string());
                        }
                    });

                // todo: This effect on water volume added?
                if let ReagentPrep::Volume(v) = &mut reagent.prep {
                    ui.label("reagant Molarity (mM):");
                    let mut val = ((*v * 1_000.) as u32).to_string();
                    if ui
                        .add(TextEdit::singleline(&mut val).desired_width(40.))
                        .changed
                    {
                        let molarity_int: u32 = val.parse().unwrap_or_default();
                        *v = molarity_int as f32 / 1_000.;
                        reagent.calc_amount(solution.total_volume);
                    }
                }

                if type_prev != reagent.type_ || prep_prev != reagent.prep {
                    reagent.calc_amount(solution.total_volume);
                }

                ui.add_space(COL_SPACING / 2.);
                ui.label("Molarity (mM):");

                // Convert to a mM integer.
                let mut val = ((reagent.molarity * 1_000.) as u32).to_string();
                if ui
                    .add(TextEdit::singleline(&mut val).desired_width(40.))
                    .changed
                {
                    let molarity_int: u32 = val.parse().unwrap_or_default();
                    reagent.molarity = molarity_int as f32 / 1_000.;
                    reagent.calc_amount(solution.total_volume);
                }

                ui.add_space(COL_SPACING / 2.);

                let result_label = if let ReagentPrep::Volume(_) = reagent.prep {
                    "Volume"
                } else {
                    "Mass"
                };

                ui.label(result_label);

                ui.add_sized(
                    [80.0, 20.0],
                    egui::Label::new(
                        RichText::new(reagent.amount_calc.to_string()).color(Color32::LIGHT_BLUE),
                    ),
                );

                ui.add_space(COL_SPACING);
                if ui.button(RichText::new("🗑").color(Color32::RED)).clicked() {
                    reagent_removed = Some(j);
                }
            });
        }

        if let Some(rem_i) = reagent_removed {
            solution.reagents.remove(rem_i);
        }

        ui.add_space(ROW_SPACING * 2.);
    }

    if let Some(rem_i) = sol_removed {
        portions.solutions.remove(rem_i);
    }
}

fn media_disp(portions: &mut PortionsState, ui: &mut Ui) {
    let type_prev = portions.media_input.clone();
    let mut run_calc = false;

    ui.horizontal(|ui| {
        ComboBox::from_id_source(3_000)
            .width(110.)
            .selected_text(portions.media_input.to_string())
            .show_ui(ui, |ui| {
                for type_ in [
                    MediaPrepInput::Plates((PlateSize::D90, 6)),
                    MediaPrepInput::Liquid(0.),
                ] {
                    ui.selectable_value(
                        &mut portions.media_input,
                        type_.clone(),
                        type_.to_string(),
                    );
                }
            });

        ui.add_space(COL_SPACING);

        if portions.media_input != type_prev {
            run_calc = true;
        }

        match &mut portions.media_input {
            MediaPrepInput::Plates((plate_size, num)) => {
                ui.label("Plate diameter:");

                let dia_prev = plate_size.clone();
                ComboBox::from_id_source(3_001)
                    .width(70.)
                    .selected_text(plate_size.to_string())
                    .show_ui(ui, |ui| {
                        for size in [
                            PlateSize::D60,
                            PlateSize::D90,
                            PlateSize::D100,
                            PlateSize::D150,
                        ] {
                            ui.selectable_value(plate_size, size, size.to_string());
                        }
                    });

                if *plate_size != dia_prev {
                    run_calc = true;
                }

                ui.label("Num plates:");

                let mut val = num.to_string();
                if ui
                    .add(TextEdit::singleline(&mut val).desired_width(40.))
                    .changed()
                {
                    *num = val.parse().unwrap_or_default();
                    run_calc = true;
                }
            }
            MediaPrepInput::Liquid(volume) => {
                ui.label("Volume (mL):");

                let mut val_ml = (*volume * 1_000.).to_string();
                if ui
                    .add(TextEdit::singleline(&mut val_ml).desired_width(50.))
                    .changed()
                {
                    let val_int: u32 = val_ml.parse().unwrap_or_default();
                    *volume = val_int as f32 / 1_000.;
                    run_calc = true;
                }
            }
        }
    });

    ui.add_space(ROW_SPACING);

    // todo: Adjust units etc A/R for higherh values.
    let result = &portions.media_result;
    ui.horizontal(|ui| {
        ui.label("Water: ");
        ui.label(format!("{:.1} mL", result.water * 1_000.));
        ui.add_space(COL_SPACING);

        ui.label("LB: ");
        ui.label(format!("{:.2} g", result.food));
        ui.add_space(COL_SPACING);

        if result.agar > 0. {
            ui.label("Agar: ");
            ui.label(format!("{:.2} g", result.agar));
            ui.add_space(COL_SPACING);
        }

        ui.label("Antibiotic (1000×): ");
        ui.label(format!("{:.1} μL", result.antibiotic * 1_000.));
    });

    // todo: Run calc immediately.

    if run_calc {
        portions.media_result = media_prep(&portions.media_input);
    }
}

pub fn portions_page(portions: &mut PortionsState, ui: &mut Ui) {
    ui.add_space(ROW_SPACING / 2.);

    ui.horizontal(|ui| {
        ui.heading("Mixing portions");
        ui.add_space(COL_SPACING);

        if ui
            .button(RichText::new("➕ Add solution").color(Color32::GOLD))
            .clicked()
        {
            portions.solutions.push(Solution {
                total_volume: DEFAULT_TOTAL_VOLUME,
                reagents: vec![Reagent::default()],
                ..Default::default()
            });
        }
    });

    ui.add_space(ROW_SPACING);

    solutions_disp(portions, ui);
    ui.add_space(ROW_SPACING);

    ui.heading("Growth media");
    media_disp(portions, ui);
}
