use eframe::egui::{Color32, RichText, TextEdit, Ui};
use egui_extras::{Column, TableBuilder};

use crate::{
    gui::{NT_CHARS, ROW_SPACING},
    State,
};

/// EGUI component for the Primer page
pub fn primer_page(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        ui.heading(RichText::new("Primer QC").color(Color32::WHITE));

        if ui.button("➕ Add primer").clicked() {
            state.ui.primer_cols.push(Default::default())
        }

        if ui.button("Tune all").clicked() {
            for data in &mut state.ui.primer_cols {
                // todo
            }
        }
    });

    ui.label("Tuning instructions: Include more of the target sequence than required on the end[s] that can be tuned. These are the\
     ends that do not define your inert, gene of interest, insertion point etc. Mark that end as tunable using the \"Tune\" button.");

    ui.add_space(ROW_SPACING);

    TableBuilder::new(ui)
        .column(Column::initial(400.).resizable(true))
        .column(Column::initial(120.).resizable(true))
        .column(Column::auto().resizable(true))
        .column(Column::auto().resizable(true))
        .column(Column::initial(40.).resizable(true))
        .column(Column::auto().resizable(true))
        .column(Column::auto().resizable(true))
        .column(Column::auto().resizable(true))
        .column(Column::remainder())
        .header(20.0, |mut header| {
            header.col(|ui| {
                ui.heading(RichText::new("Sequence (5'-> 3')").color(Color32::WHITE));
            });
            header.col(|ui| {
                ui.heading("Description");
            });
            header.col(|ui| {
                ui.heading("Len");
            });
            header.col(|ui| {
                ui.heading("Q");
            });
            header.col(|ui| {
                ui.heading("TM");
            });
            header.col(|ui| {
                ui.heading("GC %");
            });
            header.col(|ui| {
                ui.heading("3' stab");
            });
            header.col(|ui| {
                ui.heading("Cplx");
            });
            header.col(|ui| {
                ui.heading("Dimer");
            });
        })
        .body(|mut body| {
            for data in &mut state.ui.primer_cols {
                body.row(30.0, |mut row| {
                    row.col(|ui| {
                        // gui.label(make_seq_str(&col.sequence));
                        // gui.label(&col.sequence_input);
                        // let mut val = col.sequence_input;

                        ui.horizontal(|ui| {
                            if ui
                                .button(RichText::new("Tunable").color(if data.tunable_5p {
                                    Color32::GREEN
                                } else {
                                    Color32::LIGHT_GRAY
                                }))
                                .clicked()
                            {
                                data.tunable_5p = !data.tunable_5p;
                            }

                            let response = ui.add(
                                TextEdit::singleline(&mut data.sequence_input).desired_width(400.),
                            );

                            if response.changed() {
                                data.sequence_input = data
                                    .sequence_input
                                    .chars()
                                    .filter(|&c| NT_CHARS.contains(&c))
                                    .collect();
                                data.sequence_input = data.sequence_input.to_lowercase();
                                data.run_calcs();
                            }

                            if ui
                                .button(RichText::new("Tunable").color(if data.tunable_3p {
                                    Color32::GREEN
                                } else {
                                    Color32::LIGHT_GRAY
                                }))
                                .clicked()
                            {
                                data.tunable_3p = !data.tunable_3p;
                            }
                        });
                    });
                    row.col(|ui| {
                        ui.add(TextEdit::singleline(&mut data.description));
                    });

                    row.col(|ui| {
                        ui.label(data.primer.sequence.len().to_string());
                    });

                    // todo: Color-code these

                    row.col(|ui| {
                        ui.label("-");
                    });

                    row.col(|ui| {
                        let text = match &data.metrics {
                            Some(m) => &format!("{:.1}°C", m.melting_temp),
                            None => "-",
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match &data.metrics {
                            // todo: Cache this calc?
                            Some(m) => &format!("{:.0}%", m.gc_portion * 100.),
                            None => "-",
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match &data.metrics {
                            Some(m) => &format!("{}", m.gc_3p_count), // todo
                            None => "-",
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match &data.metrics {
                            Some(m) => &format!("{}", m.complexity),
                            None => "-",
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match &data.metrics {
                            Some(m) => &format!("{}", m.self_end_dimer),
                            None => "-",
                        };
                        ui.label(text);
                    });
                });
            }
        });
}
