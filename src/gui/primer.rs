use eframe::egui::{Align, Color32, Layout, RichText, TextEdit, Ui};
use egui_extras::{Column, TableBuilder};

// todo: monospace font for all seqs.
use crate::{
    gui::{page_primers_selector, PagePrimer, COL_SPACING, ROW_SPACING},
    primer::{design_amplification_primers, design_slic_fc_primers},
    util,
    util::{make_seq_str, save, seq_from_str},
    State,
};
use crate::{
    gui::{page_seq_selector, seq_view::sequence_vis, PageSeq},
    primer::{
        PrimerData,
        PrimerDirection::{self, Forward, Reverse},
        TuneSetting,
    },
    util::{get_row_ranges, seq_complement, seq_i_to_pixel},
};

const COLOR_GOOD: Color32 = Color32::GREEN;
const COLOR_MARGINAL: Color32 = Color32::GOLD;
const COLOR_BAD: Color32 = Color32::LIGHT_RED;

const DEFAULT_TRIM_AMT: usize = 32 - 20;

const TABLE_ROW_HEIGHT : f32 = 60.;

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

/// Shows below each primer sequence. Data and controls on trimming primer size for optimization.
/// Returns wheather a button was clicked.
fn primer_tune_display(data: &mut PrimerData, ui: &mut Ui) -> bool {
    // This avoids a double-mutable error
    let mut tuned = false;

    // Section for tuning primer length.
    ui.horizontal(|ui| {
        // This layout allows even spacing.
        // ui.allocate_ui(egui::Vec2::new(ui.available_width(), 0.0), |ui| {
        //     ui.with_layout(Layout::left_to_right(Align::Center), |ui| {

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

            // ui.label(&format!("({i})"));
        }

        // This section shows the trimmed sequence, with the removed parts visible to the left and right.
        ui.with_layout(Layout::left_to_right(Align::Center), |ui| {
            ui.label(RichText::new(&data.seq_removed_5p).color(Color32::GRAY));
            ui.add_space(COL_SPACING / 2.);

            if data.tunable_5p != TuneSetting::Disabled || data.tunable_3p != TuneSetting::Disabled
            {
                ui.label(
                    RichText::new(make_seq_str(&data.primer.sequence)).color(Color32::LIGHT_BLUE),
                );
            }

            ui.add_space(COL_SPACING / 2.);
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

                    // todo: We still have a crash her.e
                    if *i + 1 < data.sequence_input.len() - t5p_len {
                        *i += 1;
                    }
                    tuned = true;
                };
                // ui.label(&format!("({i})"));
            }
        });

        if tuned {
            data.run_calcs();
        }
    });
    tuned
}

fn amplification(state: &mut State, ui: &mut Ui) {
    ui.heading("Amplification");

    ui.add_space(ROW_SPACING);

    ui.label("Amplicon:");
    let response =
        ui.add(TextEdit::multiline(&mut state.ui.seq_amplicon_input).desired_width(800.));
    if response.changed() {
        state.seq_amplicon = seq_from_str(&state.ui.seq_amplicon_input);
        state.ui.seq_amplicon_input = make_seq_str(&state.seq_amplicon);
        state.sync_primer_matches(None);
    }
    ui.label(&format!("len: {}", state.ui.seq_amplicon_input.len()));

    ui.add_space(ROW_SPACING);

    ui.horizontal(|ui| {
        if ui.button("➕ Make primers").clicked() {
            // state.sync_seqs();

            if let Some(primers) = design_amplification_primers(&state.seq_amplicon) {
                let sequence_input = make_seq_str(&primers.fwd.sequence);

                let mut primer_fwd = PrimerData {
                    primer: primers.fwd,
                    sequence_input,
                    description: "Amplification Fwd".to_owned(),
                    tunable_5p: TuneSetting::Disabled,
                    tunable_3p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
                    ..Default::default()
                };

                let sequence_input = make_seq_str(&primers.rev.sequence);
                let mut primer_rev = PrimerData {
                    primer: primers.rev,
                    sequence_input,
                    description: "Amplification Rev".to_owned(),
                    tunable_5p: TuneSetting::Disabled,
                    tunable_3p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
                    ..Default::default()
                };

                primer_fwd.run_calcs();
                primer_rev.run_calcs();

                state.primer_data.extend([primer_fwd, primer_rev]);

                state.sync_primer_matches(None); // note: Not requried to run on all primers.
            }
        }
    });
}

fn primer_creation_slic_fc(state: &mut State, ui: &mut Ui) {
    ui.heading("SLIC and FastCloning");

    ui.add_space(ROW_SPACING);

    ui.label("Insert:");
    let response = ui.add(
        TextEdit::multiline(&mut state.ui.seq_insert_input).desired_width(ui.available_width()),
    );
    if response.changed() {
        state.seq_insert = seq_from_str(&state.ui.seq_insert_input);
        state.ui.seq_insert_input = make_seq_str(&state.seq_insert);
        state.sync_cloning_product();
    }
    ui.label(&format!("len: {}", state.ui.seq_insert_input.len()));

    ui.add_space(ROW_SPACING);

    ui.label("Vector:");
    let response = ui.add(
        TextEdit::multiline(&mut state.ui.seq_vector_input).desired_width(ui.available_width()),
    );
    if response.changed() {
        state.seq_vector = seq_from_str(&state.ui.seq_vector_input);
        state.ui.seq_vector_input = make_seq_str(&state.seq_vector);
        state.sync_cloning_product();
    }
    ui.label(&format!("len: {}", state.ui.seq_vector_input.len()));

    ui.horizontal(|ui| {
        let mut entry = state.insert_loc.to_string();
        let response = ui.add(TextEdit::singleline(&mut entry).desired_width(40.));
        if response.changed() {
            state.insert_loc = entry.parse().unwrap_or(0);
            state.sync_cloning_product();
        }

        ui.add_space(COL_SPACING);

        if ui.button("➕ Make cloning primers").clicked() {
            // state.sync_seqs();

            if let Some(primers) =
                design_slic_fc_primers(&state.seq_vector, &state.seq_insert, state.insert_loc)
            {
                let sequence_input = make_seq_str(&primers.insert_fwd.sequence);

                let mut insert_fwd = PrimerData {
                    primer: primers.insert_fwd,
                    sequence_input,
                    description: "SLIC Insert Fwd".to_owned(),
                    // Both ends are  tunable, since this glues the insert to the vector
                    tunable_5p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
                    tunable_3p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
                    ..Default::default()
                };

                let sequence_input = make_seq_str(&primers.insert_rev.sequence);
                let mut insert_rev = PrimerData {
                    primer: primers.insert_rev,
                    sequence_input,
                    description: "SLIC Insert Rev".to_owned(),
                    // Both ends are tunable, since this glues the insert to the vector
                    tunable_5p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
                    tunable_3p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
                    ..Default::default()
                };

                let sequence_input = make_seq_str(&primers.vector_fwd.sequence);
                let mut vector_fwd = PrimerData {
                    primer: primers.vector_fwd,
                    sequence_input,
                    description: "SLIC Vector Fwd".to_owned(),
                    // 5' is non-tunable: This is the insert location.
                    tunable_5p: TuneSetting::Disabled,
                    tunable_3p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
                    ..Default::default()
                };

                let sequence_input = make_seq_str(&primers.vector_rev.sequence);
                let mut vector_rev = PrimerData {
                    primer: primers.vector_rev,
                    sequence_input,
                    description: "SLIC Vector Rev".to_owned(),
                    tunable_5p: TuneSetting::Disabled,
                    // 3' is non-tunable: This is the insert location.
                    tunable_3p: TuneSetting::Enabled(crate::gui::primer::DEFAULT_TRIM_AMT),
                    ..Default::default()
                };
                insert_fwd.run_calcs();
                insert_rev.run_calcs();
                vector_fwd.run_calcs();
                vector_rev.run_calcs();

                state
                    .primer_data
                    .extend([insert_fwd, insert_rev, vector_fwd, vector_rev]);

                state.sync_primer_matches(None); // note: Not requried to run on all primers.
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
     ends that do not define your insert, gene of interest, insertion point etc. Mark that end as tunable using the \"Tune\" button.\
     ");

    ui.add_space(ROW_SPACING);

    if let Some(sel_i) = state.ui.primer_selected {
        ui.horizontal(|ui| {
            if sel_i + 1 > state.primer_data.len() {
                // This currently happens if deleting the final primer.
                eprintln!("Error: Exceeded primer selection len");
                state.ui.primer_selected = None;
                return;
            }

            ui.heading(&format!(
                "Selected: {}",
                &state.primer_data[sel_i].description
            ));

            ui.add_space(COL_SPACING);

            if ui.button(RichText::new("Up")).clicked() {
                // todo: Arrow icons
                if sel_i != 0 {
                    state.primer_data.swap(sel_i, sel_i - 1);
                    state.ui.primer_selected = Some(sel_i - 1);
                }
            }
            if ui.button(RichText::new("Dn")).clicked() && sel_i != state.primer_data.len() - 1 {
                state.primer_data.swap(sel_i, sel_i + 1);
                state.ui.primer_selected = Some(sel_i + 1);
            }

            if ui
                .button(RichText::new("Delete 🗑").color(Color32::RED))
                .clicked()
            {
                state.primer_data.remove(sel_i);
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
                ui.heading("Sequence (5' ⏵ 3')");
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
                ui.heading("TM").on_hover_text("Primer melting temperature, in °C. See the readme for calculations and assumptions.");
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
            header.col(|ui| {});
        })
        .body(|mut body| {
            for (i, data) in state.primer_data.iter_mut().enumerate() {
                body.row(TABLE_ROW_HEIGHT, |mut row| {
                    row.col(|ui| {
                        ui.horizontal(|ui| {
                            if ui
                                .button(RichText::new("Tun").color(if let TuneSetting::Enabled(_) = data.tunable_5p {
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
                                run_match_sync = Some(i);
                            }

                            let response = ui.add(
                                TextEdit::singleline(&mut data.sequence_input).desired_width(400.),
                            );

                            if response.changed() {
                                data.sequence_input =
                                    make_seq_str(&seq_from_str(&data.sequence_input));
                                data.run_calcs();
                                run_match_sync = Some(i);
                            }

                            if ui
                                .button(RichText::new("Tun").color(if let TuneSetting::Enabled(_) = data.tunable_3p {
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
                                run_match_sync = Some(i);
                            }
                        });

                        let updated_seq = primer_tune_display(data, ui);
                        if updated_seq {
                            run_match_sync = Some(i);
                        }
                    });

                    row.col(|ui| {
                        ui.add(TextEdit::singleline(&mut data.description));
                    });

                    row.col(|ui| {
                        ui.label(data.primer.sequence.len().to_string());
                    });

                    row.col(|ui| {
                        // todo: Cache this?
                        let num_matches = data.matches_amplification_seq.len() + data.matches_vector_with_insert.len();
                        ui.label(num_matches.to_string());
                    });

                    row.col(|ui| {
                        let text = match &data.metrics {
                            // todo: PRe-compute the * 100?
                            Some(m) => RichText::new(format!("{:.0}", m.quality_score * 100.))
                                .color(color_from_score(m.quality_score)),

                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match &data.metrics {
                            Some(m) => RichText::new(format!("{:.1}°C", m.melting_temp))
                                .color(color_from_score(m.tm_score)),

                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match &data.metrics {
                            // todo: Cache this calc?
                            Some(m) => RichText::new(format!("{:.0}%", m.gc_portion * 100.))
                                .color(color_from_score(m.gc_score)),
                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match &data.metrics {
                            Some(m) => RichText::new(format!("{}", m.gc_3p_count))
                                .color(color_from_score(m.gc_3p_score)),
                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

                    // row.col(|ui| {
                    //     let text = match &data.metrics {
                    //         Some(m) => RichText::new(format!("{}", m.complexity))
                    //             .color(color_from_score(m.complexity_score)),
                    //         None => RichText::new("-"),
                    //     };
                    //     ui.label(text);
                    // });

                    row.col(|ui| {
                        let text = match &data.metrics {
                            Some(m) => RichText::new(format!("{}", m.self_end_dimer))
                                .color(color_from_score(m.dimer_score)),
                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

                    row.col(|ui| {
                        let text = match &data.metrics {
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
        state.sync_primer_matches(run_match_sync);
    }

    ui.add_space(ROW_SPACING * 1.);

    ui.horizontal(|ui| {
        page_primers_selector(state, ui);
        ui.add_space(2. * COL_SPACING);
        page_seq_selector(state, ui);
    });

    ui.add_space(ROW_SPACING);

    match state.ui.page_seq {
        PageSeq::Edit => match state.ui.page_primer {
            PagePrimer::Amplification => {
                amplification(state, ui);
            }
            PagePrimer::SlicFc => {
                primer_creation_slic_fc(state, ui);
            }
        },
        PageSeq::View => {
            match state.ui.page_primer {
                PagePrimer::SlicFc => {

                    ui.horizontal(|ui| {
                        // todo: DRY with above
                        ui.label("Insert location: ");
                        let mut entry = state.insert_loc.to_string();
                        if ui.button("⏴").clicked() {
                            if state.insert_loc > 0 {
                                state.insert_loc -= 1;
                            }
                            state.sync_cloning_product();
                        };

                        let response = ui.add(TextEdit::singleline(&mut entry).desired_width(40.));
                        if response.changed() {
                            state.insert_loc = entry.parse().unwrap_or(0);
                            state.sync_cloning_product();
                        }

                        if ui.button("⏵").clicked() {
                            if state.insert_loc + 1 < state.seq_vector.len() {
                                state.insert_loc += 1;
                            }
                            state.sync_cloning_product();
                        };
                    });
                }

                _ => (),
            }

            sequence_vis(&state, ui);
        }
    }
}
