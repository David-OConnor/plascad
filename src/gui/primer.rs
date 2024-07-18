use eframe::egui::{Align, Align2, Color32, FontFamily, FontId, Frame, Layout, lerp, Pos2, pos2, Rect, remap, RichText, Shape, TextEdit, Ui, vec2};
use eframe::emath::RectTransform;
use eframe::epaint::{Fonts, PathStroke};
use egui_extras::{Column, TableBuilder};

use crate::{
    gui::{COL_SPACING, page_primers_selector, PagePrimerCreation, ROW_SPACING},
    primer::{design_amplification_primers, design_slic_fc_primers},
    State,
    util::{make_seq_str, save, seq_from_str},
};
use crate::primer::{PrimerData, TuneSetting};

const COLOR_GOOD: Color32 = Color32::GREEN;
const COLOR_MARGINAL: Color32 = Color32::GOLD;
const COLOR_BAD: Color32 = Color32::LIGHT_RED;

const DEFAULT_TRIM_AMT: usize = 32 - 20;

// const TM_IDEAL: f32 = 59.; // todo: Fill thi sin
//
// const THRESHOLDS_TM: (f32, f32) = (59., 60.);
// const THRESHOLDS_GC: (f32, f32) = (59., 60.);

// todo: Move this graphics drawing code to a new module, A/R.
/// Draw the sequence with primers, insertion points, and other data visible, A/R
fn sequence_vis(state: &State, ui: &mut Ui) {
    Frame::canvas(ui.style()).show(ui, |ui| {
        let mut shapes = vec![];

        let color = Color32::BLUE;

        let time = ui.input(|i| i.time);
        let speed = 1.5;
        let mode = 2.;
        let n = 120;

        let desired_size = ui.available_width() * vec2(1.0, 0.15);
        let (_id, rect) = ui.allocate_space(desired_size);

        let ctx = ui.ctx();

        let to_screen =
            RectTransform::from_to(Rect::from_x_y_ranges(0.0..=1.0, -1.0..=1.0), rect);

        let points: Vec<Pos2> = (0..=n)
            .map(|i| {
                let t = i as f64 / (n as f64);
                let amp = (time * speed * mode).sin() / mode;
                let y = amp * (t * std::f64::consts::TAU / 2.0 * mode).sin();
                to_screen * pos2(t as f32, y as f32)
            })
            .collect();

        let thickness = 10.0 / mode as f32;
        // shapes.push(Shape::line(
        //     points,
        //     // if self.colors {
        //     //     PathStroke::new_uv(thickness, move |rect, p| {
        //     //         let t = remap(p.x, rect.x_range(), -1.0..=1.0).abs();
        //     //         let center_color = hex_color!("#5BCEFA");
        //     //         let outer_color = hex_color!("#F5A9B8");
        //     //
        //     //         Color32::from_rgb(
        //     //             lerp(center_color.r() as f32..=outer_color.r() as f32, t) as u8,
        //     //             lerp(center_color.g() as f32..=outer_color.g() as f32, t) as u8,
        //     //             lerp(center_color.b() as f32..=outer_color.b() as f32, t) as u8,
        //     //         )
        //     //     })
        //     // } else {
        //     PathStroke::new(thickness, color)
        //     // },
        // ));

        // shapes.push(Shape::Text("actga"));
        let seq_text = ctx.fonts(|fonts| {
            Shape::text(
                fonts,
                Pos2::new(100., 600.),
                Align2::CENTER_CENTER,
                "actgggaacc",
                FontId::new(20., FontFamily::Proportional),
                Color32::LIGHT_BLUE,
                // ctx.style().visuals.text_color(),
            )
        });

        shapes.push(seq_text);

        ui.painter().extend(shapes);
    });

    ui.add_space(ROW_SPACING);
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

                    // todo: We still have a crash her.e
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
    ui.label(&format!("len: {}", state.ui.seq_amplicon_input.len()));

    ui.add_space(ROW_SPACING);

    ui.horizontal(|ui| {
        if ui.button("➕ Make primers").clicked() {
            state.sync_seqs();

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
            }
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
    ui.label(&format!("len: {}", state.ui.seq_insert_input.len()));

    ui.add_space(ROW_SPACING);

    ui.label("Vector:");
    let response = ui.add(TextEdit::multiline(&mut state.ui.seq_vector_input).desired_width(800.));
    if response.changed() {
        state.ui.seq_vector_input = make_seq_str(&seq_from_str(&state.ui.seq_vector_input));
    }
    ui.label(&format!("len: {}", state.ui.seq_vector_input.len()));

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
                    tunable_5p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
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

    if let Some(sel_i) = state.ui.primer_selected {
        ui.horizontal(|ui| {
            ui.heading(&format!(
                "Selected: {}",
                &state.primer_data[sel_i].description
            ));

            ui.add_space(COL_SPACING);

            if ui
                .button(RichText::new("Up"))
                .clicked()
            {
                // todo: Arrow icons
                if sel_i != 0 {
                    state.primer_data.swap(sel_i, sel_i - 1);
                    state.ui.primer_selected = Some(sel_i - 1);
                }
            }
            if ui
                .button(RichText::new("Dn"))
                .clicked() && sel_i != state.primer_data.len() - 1 {
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

    TableBuilder::new(ui)
        .column(Column::initial(700.).resizable(true))
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
                ui.heading("3'GC").on_hover_text("3' end stability: The number of G or C nucleotides in the last 5 nucleotides.");
            });
            header.col(|ui| {
                ui.heading("Cplx").on_hover_text("Sequence complexity. See the readme for calculations and assumptions.");
            });
            header.col(|ui| {
                ui.heading("Dmr").on_hover_text("Potential of forming a self-end dimer. See the readme for calculations and assumptions.");
            });
            header.col(|ui| {
                ui.heading("Rep").on_hover_text("Count of repeats of a single or double nt sequence >4 in a row.");
            });

            // For selecting the row.
            header.col(|ui| {});
        })
        .body(|mut body| {
            for (i, data) in state.primer_data.iter_mut().enumerate() {
                body.row(30.0, |mut row| {
                    row.col(|ui| {
                        // gui.label(make_seq_str(&col.sequence));
                        // gui.label(&col.sequence_input);
                        // let mut val = col.sequence_input;

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

                    row.col(|ui| {
                        let text = match &data.metrics {
                            Some(m) => RichText::new(format!("{}", m.complexity))
                                .color(color_from_score(m.complexity_score)),
                            None => RichText::new("-"),
                        };
                        ui.label(text);
                    });

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

    ui.add_space(ROW_SPACING * 3.);

    // todo: Only if you have  sequence of some sort
    sequence_vis(&state, ui);

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
