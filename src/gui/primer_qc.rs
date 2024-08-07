//! This module contains code to the primer editor, QC etc.

use eframe::egui::{Align, Color32, Layout, RichText, TextEdit, Ui};
use egui_extras::{Column, TableBuilder};

use crate::{
    gui::{COL_SPACING, ROW_SPACING},
    primer::{make_amplification_primers, Primer, PrimerData, TuneSetting},
    sequence::{seq_from_str, seq_to_str},
    IonConcentrations, State,
};

const TABLE_ROW_HEIGHT: f32 = 60.;

const COLOR_GOOD: Color32 = Color32::GREEN;
const COLOR_MARGINAL: Color32 = Color32::GOLD;
const COLOR_BAD: Color32 = Color32::LIGHT_RED;

pub const DEFAULT_TRIM_AMT: usize = 32 - 20;

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

/// Allows editing ion concentration, including float manip. Return if the response changed,
/// so we can redo TM calcs downstream.
fn ion_edit(val: &mut f32, label: &str, ui: &mut Ui) -> bool {
    ui.label(label);

    let mut v = format!("{:.1}", val);
    let response = ui.add(TextEdit::singleline(&mut v).desired_width(30.));

    if response.changed() {
        *val = v.parse().unwrap_or(0.);
        true
    } else {
        false
    }
}

pub fn primer_details(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        let add_btn = ui
            .button("➕ Add primer")
            .on_hover_text("Adds a primer to the list below. Ctrl + A");
        if add_btn.clicked() {
            state.generic.primers.push(Default::default())
        }

        if ui
            .button("➕ Make whole seq primers")
            .on_hover_text("Adds a primer pair that amplify the entire loaded sequence.")
            .clicked()
        {
            make_amplification_primers(state);
        }

        let mut sync_primer_matches = false; // Prevents a double-borrow error.
        if ui.button("Tune all").clicked() {
            for primer in &mut state.generic.primers {
                primer.tune(&state.ion_concentrations);
                sync_primer_matches = true;
            }
        }

        if sync_primer_matches {
            state.sync_primer_matches(None);
        }

        ui.add_space(COL_SPACING * 2.);

        ui.add_space(2. * COL_SPACING);

        ui.heading("Ions: (mMol)");

        if ion_edit(&mut state.ion_concentrations.monovalent, "Na+ and K+", ui)
            || ion_edit(&mut state.ion_concentrations.divalent, "mg2+", ui)
            || ion_edit(&mut state.ion_concentrations.dntp, "dNTP", ui)
            || ion_edit(&mut state.ion_concentrations.primer, "primer (nM)", ui)
        {
            for primer in &mut state.generic.primers {
                primer.run_calcs(&state.ion_concentrations); // Note: We only need to run the TM calc.
            }
        }

        // if ui.button("Load").clicked() {}

        // // todo: Temp. Find a better way.
        // if ui.button("Sync primer disp").clicked() {
        //     for p_data in &mut state.primer_data {
        //         p_data.matches_amplification_seq = p_data.primer.match_to_seq(&state.seq_amplicon);
        //         p_data.matches_slic_insert = p_data.primer.match_to_seq(&state.seq_insert);
        //         p_data.matches_slic_vector = p_data.primer.match_to_seq(&state.seq_vector);
        //     }
        // }
    });

    ui.label("Tuning instructions: Include more of the target sequence than required on the end[s] that can be tuned. These are the \
     ends that do not define your insert, gene of interest, insertion point etc. Mark that end as tunable using the \"T\" button. \
To learn about a table column, mouse over it.");

    ui.add_space(ROW_SPACING);

    if let Some(sel_i) = state.ui.primer_selected {
        ui.horizontal(|ui| {
            if sel_i + 1 > state.generic.primers.len() {
                // This currently happens if deleting the bottom-most primer.
                // If so, select the primer above it.
                state.ui.primer_selected = if !state.generic.primers.is_empty() {
                    Some(state.generic.primers.len() - 1)
                } else {
                    None
                };
                return;
            }

            ui.heading(&format!(
                "Selected: {}",
                &state.generic.primers[sel_i].description
            ));

            ui.add_space(COL_SPACING);

            if ui.button(RichText::new("Up")).clicked() {
                // todo: Arrow icons
                if sel_i != 0 {
                    state.generic.primers.swap(sel_i, sel_i - 1);
                    state.ui.primer_selected = Some(sel_i - 1);
                }
            }
            if ui.button(RichText::new("Dn")).clicked() && sel_i != state.generic.primers.len() - 1
            {
                state.generic.primers.swap(sel_i, sel_i + 1);
                state.ui.primer_selected = Some(sel_i + 1);
            }

            if ui
                .button(RichText::new("Delete 🗑").color(Color32::RED))
                .clicked()
            {
                state.generic.primers.remove(sel_i);
            }

            if ui
                .button(RichText::new("Deselect").color(Color32::GOLD))
                .clicked()
            {
                state.ui.primer_selected = None;
            }
        });

        ui.add_space(ROW_SPACING);
    }

    let mut run_match_sync = None; // Avoids a double-mutation error.

    TableBuilder::new(ui)
        .column(Column::initial(600.).resizable(true)) // Sequence
        .column(Column::initial(160.).resizable(true)) // Description
        .column(Column::auto().resizable(true))// Len
        .column(Column::auto().resizable(true))// Matches
        .column(Column::auto().resizable(true))// Quality
        .column(Column::initial(40.).resizable(true)) // TM
        .column(Column::initial(36.).resizable(true))// GC %
        .column(Column::auto().resizable(true))// 3' GC content
        .column(Column::auto().resizable(true))// Complexity
        .column(Column::auto().resizable(true))// Dimer formation
        .column(Column::auto().resizable(true))  // Repeats
        .column(Column::remainder())
        .header(20.0, |mut header| {
            header.col(|ui| {
                ui.heading("Primer sequence (5' ⏵ 3')");
            });
            header.col(|ui| {
                ui.heading("Description");
            });
            header.col(|ui| {
                ui.heading("Len").on_hover_text("Number of nucleotides in the (tuned, if applicable) primer");
            });
            header.col(|ui| {
                ui.heading("Mt").on_hover_text("Number of matches with the target sequence");
            });
            header.col(|ui| {
                ui.heading("Qual").on_hover_text("Overall primer quality. This is an abstract estimate, taking all other listed factors into account.");
            });
            header.col(|ui| {
                ui.heading("TM").on_hover_text("Primer melting temperature, in °C. Calculated using base a base stacking method, where\
                 enthalpy and entropy of neighboring base pairs are added. See the readme for calculations and assumptions.");
            });
            header.col(|ui| {
                ui.heading("GC").on_hover_text("The percentage of nucleotides that are C or G.");
            });
            header.col(|ui| {
                ui.heading("3'GC").on_hover_text("3' end stability: The number of G or C nucleotides in the last 5 nucleotides. 2-3 is ideal. Sources differ on if 4 is acceptable.");
            });
            // header.col(|ui| {
            //     ui.heading("Cplx").on_hover_text("Sequence complexity. See the readme for calculations and assumptions.");
            // });
            header.col(|ui| {
                ui.heading("Dmr").on_hover_text("Potential of forming a self-end dimer. See the readme for calculations and assumptions.");
            });
            header.col(|ui| {
                ui.heading("Rep").on_hover_text("Count of repeats of a single or double nt sequence >4 in a row, and count of triplet \
                repeats anywhere in the sequence.");
            });

            // For selecting the row.
            header.col(|_ui| {});
        })
        .body(|mut body| {
            for (i, primer) in state.generic.primers.iter_mut().enumerate() {
                body.row(TABLE_ROW_HEIGHT, |mut row| {
                    row.col(|ui| {
                        ui.horizontal(|ui| {
                            if ui
                                .button(RichText::new("T").color(if let TuneSetting::Enabled(_) = primer.volatile.tunable_5p {
                                    Color32::GREEN
                                } else {
                                    Color32::LIGHT_GRAY
                                }))
                                .clicked()
                            {
                                primer.volatile.tunable_5p.toggle();
                                if primer.volatile.tunable_5p == TuneSetting::Disabled {
                                    primer.run_calcs(&state.ion_concentrations); // To re-sync the sequence without parts removed.
                                }
                                run_match_sync = Some(i);
                            }

                            let response = ui.add(
                                TextEdit::singleline(&mut primer.volatile.sequence_input).desired_width(400.),
                            );

                            if response.changed() {
                                primer.volatile.sequence_input =
                                    seq_to_str(&seq_from_str(&primer.volatile.sequence_input));
                                primer.run_calcs(&state.ion_concentrations);
                                run_match_sync = Some(i);
                            }

                            if ui
                                .button(RichText::new("T").color(if let TuneSetting::Enabled(_) = primer.volatile.tunable_3p {
                                    Color32::GREEN
                                } else {
                                    Color32::LIGHT_GRAY
                                }))
                                .clicked()
                            {
                                primer.volatile.tunable_3p.toggle();
                                if primer.volatile.tunable_3p == TuneSetting::Disabled {
                                    primer.run_calcs(&state.ion_concentrations); // To re-sync the sequence without parts removed.
                                }
                                run_match_sync = Some(i);
                            }

                            ui.add_space(COL_SPACING);

                            if primer.volatile.tunable_3p != TuneSetting::Disabled || primer.volatile.tunable_5p != TuneSetting::Disabled {
                                if ui
                                    .button(RichText::new("Tune")).on_hover_text("Tune selected ends for this primer").clicked()
                                {
                                    primer.tune(&state.ion_concentrations);
                                    run_match_sync = Some(i);
                                }
                            }
                        });

                        let updated_seq = primer_tune_display(primer, &state.ion_concentrations, ui);
                        if updated_seq {
                            run_match_sync = Some(i);
                        }
                    });

                    row.col(|ui| {
                        ui.add(TextEdit::singleline(&mut primer.description));
                    });

                    row.col(|ui| {
                        ui.label(primer.sequence.len().to_string());
                    });

                    row.col(|ui| {
                        // todo: Cache this?
                        // let num_matches = data.matches_seq.len() + data.matches_vector_with_insert.len();
                        let num_matches = primer.volatile.matches_seq.len() +  primer.volatile.matches_seq.len();
                        ui.label(num_matches.to_string());
                    });

                    row.col(|ui| {
                        let text = match & primer.volatile.metrics {
                            // todo: PRe-compute the * 100?
                            Some(m) => RichText::new(format!("{:.0}", m.quality_score * 100.))
                                .color(color_from_score(m.quality_score)),

                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match & primer.volatile.metrics {
                            Some(m) => RichText::new(format!("{:.1}°C", m.melting_temp))
                                .color(color_from_score(m.tm_score)),

                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match & primer.volatile.metrics {
                            // todo: Cache this calc?
                            Some(m) => RichText::new(format!("{:.0}%", m.gc_portion * 100.))
                                .color(color_from_score(m.gc_score)),
                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match & primer.volatile.metrics {
                            Some(m) => RichText::new(format!("{}", m.gc_3p_count))
                                .color(color_from_score(m.gc_3p_score)),
                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match & primer.volatile.metrics {
                            Some(m) => RichText::new(format!("{}", m.self_end_dimer))
                                .color(color_from_score(m.dimer_score)),
                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match & primer.volatile.metrics {
                            Some(m) => RichText::new(format!("{}", m.repeats))
                                .color(color_from_score(m.repeats_score)),
                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let mut selected = false;

                        if let Some(sel_i) = state.ui.primer_selected {
                            if sel_i == i {
                                selected = true
                            }
                        }

                        if selected {
                            if ui.button(RichText::new("🔘").color(Color32::GREEN)).clicked() {
                                state.ui.primer_selected = None;
                            }
                        } else if ui.button("🔘").clicked() {
                            state.ui.primer_selected = Some(i);
                        }
                    });
                });
            }
        });

    if run_match_sync.is_some() {
        state.sync_seq_related(run_match_sync);
    }
}

/// Shows below each primer sequence. Data and controls on trimming primer size for optimization.
/// Returns wheather a button was clicked.
fn primer_tune_display(
    primer: &mut Primer,
    ion_concentrations: &IonConcentrations,
    ui: &mut Ui,
) -> bool {
    // This avoids a double-mutable error
    let mut tuned = false;

    // Section for tuning primer length.
    ui.horizontal(|ui| {
        // This layout allows even spacing.
        // ui.allocate_ui(egui::Vec2::new(ui.available_width(), 0.0), |ui| {
        //     ui.with_layout(Layout::left_to_right(Align::Center), |ui| {

        if let TuneSetting::Enabled(i) = &mut primer.volatile.tunable_5p {
            ui.label("5'");
            if ui.button("⏴").clicked() {
                if *i > 0 {
                    *i -= 1;
                }
                tuned = true;
            };
            if ui.button("⏵").clicked() {
                let t3p_len = match primer.volatile.tunable_3p {
                    TuneSetting::Enabled(t) => t,
                    _ => 0,
                };
                if *i + 1 < primer.volatile.sequence_input.len() - t3p_len {
                    *i += 1;
                }
                tuned = true;
            };

            // ui.label(&format!("({i})"));
        }

        // This section shows the trimmed sequence, with the removed parts visible to the left and right.
        ui.with_layout(Layout::left_to_right(Align::Center), |ui| {
            ui.label(RichText::new(&primer.volatile.seq_removed_5p).color(Color32::GRAY));
            ui.add_space(COL_SPACING / 2.);

            if primer.volatile.tunable_5p != TuneSetting::Disabled
                || primer.volatile.tunable_3p != TuneSetting::Disabled
            {
                ui.label(RichText::new(seq_to_str(&primer.sequence)).color(Color32::LIGHT_BLUE));
            }

            ui.add_space(COL_SPACING / 2.);
            ui.label(RichText::new(&primer.volatile.seq_removed_3p).color(Color32::GRAY));
        });

        // Note: We need to reverse the item order for this method of right-justifying to work.
        // This is kind of OK with the intent here though.
        ui.with_layout(Layout::right_to_left(Align::Max), |ui| {
            if let TuneSetting::Enabled(i) = &mut primer.volatile.tunable_3p {
                ui.label("3'");

                if ui.button("⏵").clicked() {
                    if *i > 0 {
                        *i -= 1;
                    }
                    tuned = true;
                };
                if ui.button("⏴").clicked() {
                    let t5p_len = match primer.volatile.tunable_5p {
                        TuneSetting::Enabled(t) => t,
                        _ => 0,
                    };

                    // todo: We still have a crash her.e
                    if *i + 1 < primer.volatile.sequence_input.len() - t5p_len {
                        *i += 1;
                    }
                    tuned = true;
                };
                // ui.label(&format!("({i})"));
            }
        });

        if tuned {
            primer.run_calcs(ion_concentrations);
        }
    });
    tuned
}
