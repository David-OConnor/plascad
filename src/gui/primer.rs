use bincode::{Decode, Encode};
use eframe::egui::{Align, Color32, Direction, Layout, RichText, TextEdit, Ui};
use egui_extras::{Column, TableBuilder};

use crate::{
    gui::{page_primers_selector, PagePrimerCreation, COL_SPACING, ROW_SPACING},
    primer::{design_slic_fc_primers, Primer, PrimerMetrics},
    util::{make_seq_str, save, seq_from_str},
    State,
};

const COLOR_GOOD: Color32 = Color32::GREEN;
const COLOR_MARGINAL: Color32 = Color32::GOLD;
const COLOR_BAD: Color32 = Color32::LIGHT_RED;

// const TM_IDEAL: f32 = 59.; // todo: Fill thi sin
//
// const THRESHOLDS_TM: (f32, f32) = (59., 60.);
// const THRESHOLDS_GC: (f32, f32) = (59., 60.);

#[derive(Clone, Copy, PartialEq, Encode, Decode)]
enum TuneSetting {
    Disabled,
    /// Inner: Offset index; this marks the distance from the respective ends that the sequence is attentuated to.
    Enabled(usize),
}

impl Default for TuneSetting {
    fn default() -> Self {
        Self::Disabled
    }
}

impl TuneSetting {
    pub fn toggle(&mut self) {
        *self = match self {
            Self::Disabled => Self::Enabled(0),
            _ => Self::Disabled,
        }
    }
}

#[derive(Default, Encode, Decode)]
pub struct PrimerData {
    /// This primer excludes nts past the tuning offsets.
    pub primer: Primer,
    /// Editing is handled using this string; we convert the string to our nucleotide sequence as needed.
    /// This includes nts past the tuning offsets.
    pub sequence_input: String,
    pub description: String,
    pub metrics: Option<PrimerMetrics>,
    // tunable_end: TunableEnd,
    /// These fields control if a given primer end is fixed (Eg marking the start of an insert,
    /// marking the insert point in a vector etc) or if we can tune its length to optimize the primer.
    pub tunable_5p: TuneSetting,
    pub tunable_3p: TuneSetting,
    /// These seq_removed fields are redundant with primer and tune settings. We use them to cache
    /// the actual sequences that are removed for display purposes.
    // seq_removed_5p: Seq,
    // seq_removed_3p: Seq
    pub seq_removed_5p: String,
    pub seq_removed_3p: String,
}

impl PrimerData {
    /// Perform calculations on primer quality and related data. Run this when the sequence changes,
    /// the tuning values change etc.
    pub fn run_calcs(&mut self) {
        let full_len = self.sequence_input.len();
        let mut start = 0;
        let mut end = full_len;

        if let TuneSetting::Enabled(i) = self.tunable_5p {
            start = i;
        }

        if let TuneSetting::Enabled(i) = self.tunable_3p {
            end = if i > full_len {
                // Prevents an overrun.
                0
            } else {
                full_len - i
            };
        }

        if start > end || start + 1 > self.sequence_input.len() {
            start = 0;
            end = full_len
        }

        self.primer.sequence = seq_from_str(&self.sequence_input[start..end]);
        self.metrics = self.primer.calc_metrics();

        self.seq_removed_5p = self.sequence_input[..start].to_owned();
        self.seq_removed_3p = self.sequence_input[end..].to_owned();
    }
}

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

/// Shows below each primer sequence. Data and controls on trimming primer size for optimization.
fn primer_tune_display(data: &mut PrimerData, ui: &mut Ui) {
    // Section for tuning primer length.
    ui.horizontal(|ui| {
        // This layout allows even spacing.
        // ui.allocate_ui(egui::Vec2::new(ui.available_width(), 0.0), |ui| {
        //     ui.with_layout(Layout::left_to_right(Align::Center), |ui| {
        // This avoids a double-mutable error
        let mut tuned = false;

        if let TuneSetting::Enabled(i) = &mut data.tunable_5p {
            ui.label("5'");
            if ui.button("⏴").clicked() {
                if *i > 0 {
                    *i -= 1;
                }
                tuned = true;
            };
            if ui.button("⏵").clicked() {
                let t3p_len = match data.tunable_3p {
                    TuneSetting::Enabled(t) => t,
                    _ => 0,
                };
                if *i + 1 < data.sequence_input.len() - t3p_len {
                    *i += 1;
                }
                tuned = true;
            };

            ui.label(&format!("({i})"));
        }

        // This section shows the trimmed sequence, with the removed parts visible to the left and right.
        ui.with_layout(Layout::left_to_right(Align::Center), |ui| {
            ui.label(RichText::new(&data.seq_removed_5p).color(Color32::GRAY));
            ui.add_space(COL_SPACING);

            if data.tunable_5p != TuneSetting::Disabled || data.tunable_3p != TuneSetting::Disabled
            {
                ui.label(
                    RichText::new(make_seq_str(&data.primer.sequence)).color(Color32::LIGHT_BLUE),
                );
            }

            ui.add_space(COL_SPACING);
            ui.label(RichText::new(&data.seq_removed_3p).color(Color32::GRAY));
        });

        // Note: We need to reverse the item order for this method of right-justifying to work.
        // This is kind of OK with the intent here though.
        ui.with_layout(Layout::right_to_left(Align::Max), |ui| {
            if let TuneSetting::Enabled(i) = &mut data.tunable_3p {
                ui.label("3'");

                if ui.button("⏵").clicked() {
                    if *i > 0 {
                        *i -= 1;
                    }
                    tuned = true;
                };
                if ui.button("⏴").clicked() {
                    let t5p_len = match data.tunable_5p {
                        TuneSetting::Enabled(t) => t,
                        _ => 0,
                    };

                    if *i + 1 < data.sequence_input.len() - t5p_len {
                        *i += 1;
                    }
                    tuned = true;
                };
                ui.label(&format!("({i})"));
            }
        });

        if tuned {
            data.run_calcs();
        }
        // });
    });
}

fn amplification(state: &mut State, ui: &mut Ui) {
    ui.heading("Amplification");

    ui.add_space(ROW_SPACING);

    ui.label("Amplicon:");
    let response =
        ui.add(TextEdit::multiline(&mut state.ui.seq_amplicon_input).desired_width(800.));
    if response.changed() {
        state.ui.seq_amplicon_input = make_seq_str(&seq_from_str(&state.ui.seq_amplicon_input));
    }

    ui.add_space(ROW_SPACING);

    ui.horizontal(|ui| {
        if ui.button("➕ Make primers").clicked() {
            state.sync_seqs();

            // if let Some(primers) =
            //     design_amplification_primers(&state.seq_amplicon)
            // {
            //     let sequence_input = make_seq_str(&primers.fwd.sequence);
            //
            //     let mut primer_fwd = PrimerData {
            //         primer: primers.insert_fwd,
            //         sequence_input,
            //         description: "Amplification Fwd".to_owned(),
            //         // Both ends are  tunable, since this glues the insert to the vector
            //         tunable_5p: TuneSetting::Enabled(0),
            //         tunable_3p: TuneSetting::Enabled(0),
            //         ..Default::default()
            //     };
            //
            //     let sequence_input = make_seq_str(&primers.rev.sequence);
            //     let mut primer_rev = PrimerData {
            //         primer: primers.insert_rev,
            //         sequence_input,
            //         description: "Amplification Rev".to_owned(),
            //         // Both ends are tunable, since this glues the insert to the vector
            //         tunable_5p: TuneSetting::Enabled(0),
            //         tunable_3p: TuneSetting::Enabled(0),
            //         ..Default::default()
            //     };
            //
            //     primer_fwd.run_calcs();
            //     primer_rev.run_calcs();
            //
            //     state
            //         .primer_data
            //         .extend([primer_fwd, primer_rev]);
            // }
        }
    });
}

fn primer_creation_slic_fc(state: &mut State, ui: &mut Ui) {
    ui.heading("SLIC and FastCloning");

    ui.add_space(ROW_SPACING);

    ui.label("Insert:");
    let response = ui.add(TextEdit::multiline(&mut state.ui.seq_insert_input).desired_width(800.));
    if response.changed() {
        state.ui.seq_insert_input = make_seq_str(&seq_from_str(&state.ui.seq_insert_input));
    }

    ui.add_space(ROW_SPACING);

    ui.label("Vector:");
    let response = ui.add(TextEdit::multiline(&mut state.ui.seq_vector_input).desired_width(800.));
    if response.changed() {
        state.ui.seq_vector_input = make_seq_str(&seq_from_str(&state.ui.seq_vector_input));
    }

    ui.horizontal(|ui| {
        let mut entry = state.insert_loc.to_string();
        let response = ui.add(TextEdit::singleline(&mut entry).desired_width(40.));
        if response.changed() {
            state.insert_loc = entry.parse().unwrap_or(0);
        }

        ui.add_space(COL_SPACING);

        if ui.button("➕ Make cloning primers").clicked() {
            state.sync_seqs();

            if let Some(primers) =
                design_slic_fc_primers(&state.seq_vector, &state.seq_insert, state.insert_loc)
            {
                let sequence_input = make_seq_str(&primers.insert_fwd.sequence);

                let mut insert_fwd = PrimerData {
                    primer: primers.insert_fwd,
                    sequence_input,
                    description: "SLIC Insert Fwd".to_owned(),
                    // Both ends are  tunable, since this glues the insert to the vector
                    tunable_5p: TuneSetting::Enabled(0),
                    tunable_3p: TuneSetting::Enabled(0),
                    ..Default::default()
                };

                let sequence_input = make_seq_str(&primers.insert_rev.sequence);
                let mut insert_rev = PrimerData {
                    primer: primers.insert_rev,
                    sequence_input,
                    description: "SLIC Insert Rev".to_owned(),
                    // Both ends are tunable, since this glues the insert to the vector
                    tunable_5p: TuneSetting::Enabled(0),
                    tunable_3p: TuneSetting::Enabled(0),
                    ..Default::default()
                };

                let sequence_input = make_seq_str(&primers.vector_fwd.sequence);
                let mut vector_fwd = PrimerData {
                    primer: primers.vector_fwd,
                    sequence_input,
                    description: "SLIC Vector Fwd".to_owned(),
                    // 5' is non-tunable: This is the insert location.
                    tunable_5p: TuneSetting::Disabled,
                    tunable_3p: TuneSetting::Enabled(0),
                    ..Default::default()
                };

                let sequence_input = make_seq_str(&primers.vector_rev.sequence);
                let mut vector_rev = PrimerData {
                    primer: primers.vector_rev,
                    sequence_input,
                    description: "SLIC Vector Rev".to_owned(),
                    tunable_5p: TuneSetting::Enabled(0),
                    // 3' is non-tunable: This is the insert location.
                    tunable_3p: TuneSetting::Disabled,
                    ..Default::default()
                };
                insert_fwd.run_calcs();
                insert_rev.run_calcs();
                vector_fwd.run_calcs();
                vector_rev.run_calcs();

                state
                    .primer_data
                    .extend([insert_fwd, insert_rev, vector_fwd, vector_rev]);
            }
        }
    });
}

/// EGUI component for the Primer page
pub fn primer_page(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        ui.heading(RichText::new("Primer QC").color(Color32::WHITE));
        ui.add_space(COL_SPACING);

        let add_btn = ui
            .button("➕ Add primer")
            .on_hover_text("Adds a primer to the list below. Ctrl + A");
        if add_btn.clicked() {
            state.primer_data.push(Default::default())
        }

        if ui.button("Tune all").clicked() {
            for data in &mut state.primer_data {
                // todo
            }
        }

        ui.add_space(COL_SPACING * 2.);

        // todo: Ctrl + S as well.
        if ui
            .button("Save")
            .on_hover_text("Save primer data. Ctrl + S")
            .clicked()
        {
            if let Err(e) = save("plasmid_tools.save", state) {
                println!("Error saving: {e}");
            }
        }

        // if ui.button("Load").clicked() {}
    });

    ui.label("Tuning instructions: Include more of the target sequence than required on the end[s] that can be tuned. These are the \
     ends that do not define your insert, gene of interest, insertion point etc. Mark that end as tunable using the \"Tune\" button.\
     ");

    ui.add_space(ROW_SPACING);

    TableBuilder::new(ui)
        .column(Column::initial(400.).resizable(true))
        .column(Column::initial(160.).resizable(true))
        .column(Column::auto().resizable(true))
        .column(Column::auto().resizable(true))
        .column(Column::initial(40.).resizable(true))
        .column(Column::initial(36.).resizable(true))
        .column(Column::auto().resizable(true))
        .column(Column::auto().resizable(true))
        .column(Column::auto().resizable(true))
        .column(Column::auto().resizable(true))
        .column(Column::remainder())
        .header(20.0, |mut header| {
            header.col(|ui| {
                ui.heading("Sequence (5' ⏵ 3')");
            });
            header.col(|ui| {
                ui.heading("Description");
            });
            header.col(|ui| {
                ui.heading("Len").on_hover_text("Number of nucleotides in the (tuned, if applicable) primer");
            });
            header.col(|ui| {
                ui.heading("Qual").on_hover_text("Overall primer quality. This is an abstract estimate, taking all other listed factors into account.");
            });
            header.col(|ui| {
                ui.heading("TM").on_hover_text("Primer melting temperature, in °C. See the readme for calculations and assumptions.");
            });
            header.col(|ui| {
                ui.heading("GC").on_hover_text("The percentage of nucleotides that are C or G.");
            });
            header.col(|ui| {
                ui.heading("3' stab").on_hover_text("3' end stability: The number of G or C nucleotides in the last 5 nucleotides.");
            });
            header.col(|ui| {
                ui.heading("Cplx").on_hover_text("Sequence complexity. See the readme for calculations and assumptions.");
            });
            header.col(|ui| {
                ui.heading("Dimer").on_hover_text("Potential of forming a self-end dimer. See the readme for calculations and assumptions.");
            });
            header.col(|ui| {
                ui.heading("Repeat").on_hover_text("Count of repeats of a single or double nt sequence >4 in a row.");
            });

            // For deleting the row.
            header.col(|ui| {});
        })
        .body(|mut body| {
            let mut removed = None;
            for (i, data) in state.primer_data.iter_mut().enumerate() {
                body.row(30.0, |mut row| {
                    row.col(|ui| {
                        // gui.label(make_seq_str(&col.sequence));
                        // gui.label(&col.sequence_input);
                        // let mut val = col.sequence_input;

                        ui.horizontal(|ui| {
                            if ui
                                .button(RichText::new("Tunable").color(if let TuneSetting::Enabled(_) = data.tunable_5p {
                                    Color32::GREEN
                                } else {
                                    Color32::LIGHT_GRAY
                                }))
                                .clicked()
                            {
                                data.tunable_5p.toggle();
                                if data.tunable_5p == TuneSetting::Disabled {
                                    data.run_calcs(); // To re-sync the sequence without parts removed.
                                }
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
                                .button(RichText::new("Tunable").color(if let TuneSetting::Enabled(_) = data.tunable_3p {
                                    Color32::GREEN
                                } else {
                                    Color32::LIGHT_GRAY
                                }))
                                .clicked()
                            {
                                data.tunable_3p.toggle();
                                if data.tunable_3p == TuneSetting::Disabled {
                                    data.run_calcs(); // To re-sync the sequence without parts removed.
                                }
                            }
                        });

                        primer_tune_display(data, ui);
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

                    row.col(|ui| {
                        let text = match &data.metrics {
                            Some(m) => RichText::new(&format!("{}", m.repeats))
                                .color(color_from_score(m.repeats_score)),
                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        if ui.button("🗑").clicked() {
                            removed = Some(i);
                        }
                    });
                });
            }
            if let Some(rem_i) = removed {
                state.primer_data.remove(rem_i);
            }
        });

    ui.add_space(ROW_SPACING * 2.);

    page_primers_selector(state, ui);

    match state.ui.page_primer_creation {
        PagePrimerCreation::Amplification => {
            amplification(state, ui);
        }
        PagePrimerCreation::SlicFc => {
            primer_creation_slic_fc(state, ui);
        }
    }

    // todo: Visualizer here with the seq, the primers etc
}
