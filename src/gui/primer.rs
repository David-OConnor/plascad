use eframe::egui::{Color32, RichText, TextEdit, Ui};
use egui_extras::{Column, TableBuilder};

use crate::{
    gui::{COL_SPACING, ROW_SPACING},
    make_seq_str, seq_from_str, State,
};

const COLOR_GOOD: Color32 = Color32::GREEN;
const COLOR_MARGINAL: Color32 = Color32::GOLD;
const COLOR_BAD: Color32 = Color32::LIGHT_RED;

// const TM_IDEAL: f32 = 59.; // todo: Fill thi sin
//
// const THRESHOLDS_TM: (f32, f32) = (59., 60.);
// const THRESHOLDS_GC: (f32, f32) = (59., 60.);

/// Color scores in each category according to these thresholds. These scores should be on a scale
/// between 0 and 1.
fn color_from_score(score: f32) -> Color32 {
    const SCORE_COLOR_THRESH: (f32, f32) = (0.5, 0.8);

    if score > SCORE_COLOR_THRESH.1 {
        COLOR_GOOD
    } else if score > SCORE_COLOR_THRESH.0 {
        COLOR_MARGINAL
    } else {
        COLOR_BAD
    }
}

/// EGUI component for the Primer page
pub fn primer_page(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        ui.heading(RichText::new("Primer QC").color(Color32::WHITE));
        ui.add_space(COL_SPACING);

        if ui.button("➕ Add primer").clicked() {
            state.ui.primer_cols.push(Default::default())
        }

        if ui.button("Tune all").clicked() {
            for data in &mut state.ui.primer_cols {
                // todo
            }
        }

        ui.add_space(COL_SPACING * 2.);

        if ui.button("Save").clicked() {}

        if ui.button("Load").clicked() {}
    });

    ui.label("Tuning instructions: Include more of the target sequence than required on the end[s] that can be tuned. These are the\
     ends that do not define your inert, gene of interest, insertion point etc. Mark that end as tunable using the \"Tune\" button.");

    ui.add_space(ROW_SPACING);

    TableBuilder::new(ui)
        .column(Column::initial(400.).resizable(true))
        .column(Column::initial(160.).resizable(true))
        .column(Column::auto().resizable(true))
        .column(Column::auto().resizable(true))
        .column(Column::initial(40.).resizable(true))
        .column(Column::auto().resizable(true))
        .column(Column::auto().resizable(true))
        .column(Column::auto().resizable(true))
        .column(Column::remainder())
        .header(20.0, |mut header| {
            header.col(|ui| {
                ui.heading("Sequence (5'-> 3')");
            });
            header.col(|ui| {
                ui.heading("Description");
            });
            header.col(|ui| {
                ui.heading("Len");
            });
            header.col(|ui| {
                ui.heading("Qual");
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
                                data.sequence_input =
                                    make_seq_str(&seq_from_str(&data.sequence_input));
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

                    row.col(|ui| {
                        let text = match &data.metrics {
                            // todo: PRe-compute the * 100?
                            Some(m) => RichText::new(&format!("{:.0}", m.quality_score * 100.))
                                .color(color_from_score(m.quality_score)),

                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match &data.metrics {
                            Some(m) => RichText::new(&format!("{:.1}°C", m.melting_temp))
                                .color(color_from_score(m.tm_score)),

                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match &data.metrics {
                            // todo: Cache this calc?
                            Some(m) => RichText::new(&format!("{:.0}%", m.gc_portion * 100.))
                                .color(color_from_score(m.gc_score)),
                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match &data.metrics {
                            Some(m) => RichText::new(&format!("{}", m.gc_3p_count))
                                .color(color_from_score(m.gc_3p_score)),
                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match &data.metrics {
                            Some(m) => RichText::new(&format!("{}", m.complexity))
                                .color(color_from_score(m.complexity_score)),
                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match &data.metrics {
                            Some(m) => RichText::new(&format!("{}", m.self_end_dimer))
                                .color(color_from_score(m.dimer_score)),
                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });
                });
            }
        });

    ui.add_space(ROW_SPACING * 2.);

    ui.heading("SLIC and FastCloning");

    ui.add_space(ROW_SPACING);

    ui.heading("Insert:");
    let response = ui.add(TextEdit::multiline(&mut state.ui.seq_insert_input).desired_width(800.));
    if response.changed() {
        state.ui.seq_insert_input = make_seq_str(&seq_from_str(&state.ui.seq_insert_input));
    }

    ui.add_space(ROW_SPACING);

    ui.heading("Vector:");
    let response = ui.add(TextEdit::multiline(&mut state.ui.seq_vector_input).desired_width(800.));
    if response.changed() {
        state.ui.seq_vector_input = make_seq_str(&seq_from_str(&state.ui.seq_vector_input));
    }
}
