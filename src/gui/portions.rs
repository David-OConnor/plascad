//! UI page for mixing portions. (growth media, stock solutions etc)

use std::fmt::Display;

use eframe::{
    egui,
    egui::{Color32, ComboBox, RichText, TextEdit, Ui},
};

use crate::{
    gui::{COL_SPACING, ROW_SPACING},
    portions::{AmountCalculated, PortionsState, Reagent, ReagentPrep, ReagentType, Solution},
};
// todo: Make a non-gui portions module once this becomes unweildy.

// todo: Store solutions. Save to file, and mix solutions from other solutions.

pub fn portions_page(portions: &mut PortionsState, ui: &mut Ui) {
    ui.add_space(ROW_SPACING / 2.);

    ui.horizontal(|ui| {
        ui.heading("Mixing portions");
        ui.add_space(COL_SPACING);

        if ui
            .button(RichText::new("➕ Add solution").color(Color32::GOLD))
            .clicked()
        {
            portions.solutions.push(Solution::default());
        }
    });

    ui.add_space(ROW_SPACING);

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
                solution.reagents.push(Reagent {
                    type_: ReagentType::SodiumChloride,
                    prep: ReagentPrep::Mass,
                    molarity: 0.,
                    amount_calc: AmountCalculated::Mass(0.),
                });
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
                            ReagentType::SodiumChloride,
                            ReagentType::SodiumPhosphateMonobasic,
                            ReagentType::SodiumPhosphateDibasic,
                            ReagentType::SodiumPhosphateDibasicHeptahydrate,
                            ReagentType::TrisHcl,
                            ReagentType::Imidazole,
                            ReagentType::Custom(0.),
                        ] {
                            ui.selectable_value(&mut reagent.type_, type_, type_.to_string());
                        }

                        // ui.selectable_value(val, Reagent::Custom(_), Reagent::Custom(_).to_string());
                    });

                ui.add_sized(
                    [80.0, 20.0],
                    egui::Label::new(format!("{:.2} g/mol", reagent.type_.weight())),
                );

                let prep_prev = reagent.prep.clone();
                ComboBox::from_id_source(2000 + i * 100 + j)
                    .width(80.)
                    .selected_text(reagent.prep.to_string())
                    .show_ui(ui, |ui| {
                        for prep in [ReagentPrep::Mass, ReagentPrep::Volume(0.)] {
                            ui.selectable_value(&mut reagent.prep, prep, prep.to_string());
                        }

                        // ui.selectable_value(&mut component.prep, ReagentPrep::Mass, ReagentPrep::Mass.to_string());
                        // ui.selectable_value(&mut component.prep, ReagentPrep::Volume(0.), ReagentPrep::Volume(0.).to_string());

                        // ui.selectable_value(val, Reagent::Custom(_), Reagent::Custom(_).to_string());
                    });

                // todo: This effect on water volume added?
                if let ReagentPrep::Volume(v) = &mut reagent.prep {
                    ui.label("reagant Molarity (mM)");
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
                ui.label("Mass:");

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

    // ui.heading("Growth media");
}
