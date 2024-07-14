use eframe::{
    egui,
    egui::{Color32, Context, RichText},
};
use egui_extras::{Column, TableBuilder};

use crate::{
    make_seq_str,
    primer::{Primer, PrimerMetrics},
    seq_from_str, Seq, State,
};

pub const WINDOW_WIDTH: f32 = 1100.;
pub const WINDOW_HEIGHT: f32 = 1_000.;

pub const WINDOW_TITLE: &str = "Plasmid check";

pub const ROW_SPACING: f32 = 22.;
pub const COL_SPACING: f32 = 30.;

const NT_CHARS: [char; 8] = ['A', 'C', 'T', 'G', 'a', 'c', 't', 'g'];

#[derive(Clone, Copy)]
enum Page {
    Sequence,
    Enzymes,
    Primters,
    Features,
}

#[derive(Default)]
pub struct PrimerTableCol {
    primer: Primer,
    /// Editing is handled using this string; we convert the string to our nucleotide sequence as needed.
    sequence_input: String,
    description: String,
    metrics: Option<PrimerMetrics>,
}

impl PrimerTableCol {
    /// Perform calculations on primer quality and related data. Run this when the sequence changes.
    pub fn run_calcs(&mut self) {
        self.primer.sequence = seq_from_str(&self.sequence_input);
        self.metrics = self.primer.calc_metrics();

        println!("Metrics: {:?}", self.metrics);
    }
}

pub fn draw(state: &mut State, ctx: &Context) {
    egui::CentralPanel::default().show(ctx, |ui| {
        // todo: Consider adding status like tshark/Wifi connection, monitor mode enabled, GPS coords etc.

        egui::containers::ScrollArea::vertical().show(ui, |ui| {
            ui.heading(RichText::new("Primer QC").color(Color32::WHITE));

            if ui.button("➕").clicked() {
                state.ui.primer_cols.push(Default::default())
            }

            let primer_table = TableBuilder::new(ui)
                .column(Column::initial(300.).resizable(true))
                .column(Column::initial(120.).resizable(true))
                .column(Column::auto().resizable(true))
                .column(Column::auto().resizable(true))
                .column(Column::auto().resizable(true))
                .column(Column::auto().resizable(true))
                .column(Column::auto().resizable(true))
                .column(Column::remainder())
                .header(20.0, |mut header| {
                    header.col(|ui| {
                        ui.heading("Sequence");
                    });
                    header.col(|ui| {
                        ui.heading("Description");
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
                                // ui.label(make_seq_str(&col.sequence));
                                // ui.label(&col.sequence_input);
                                // let mut val = col.sequence_input;

                                let response =
                                    ui.add(egui::TextEdit::singleline(&mut data.sequence_input));

                                if response.changed() {
                                    data.sequence_input = data
                                        .sequence_input
                                        .chars()
                                        .filter(|&c| NT_CHARS.contains(&c))
                                        .collect();
                                    data.sequence_input = data.sequence_input.to_lowercase();
                                    data.run_calcs();
                                }
                            });
                            row.col(|ui| {
                                ui.add(egui::TextEdit::singleline(&mut data.description));
                            });

                            row.col(|ui| {
                                ui.label("-");
                            });

                            row.col(|ui| {
                                let text = match &data.metrics {
                                    Some(m) => &format!("{}°C", m.melting_temp),
                                    None => "-",
                                };
                                ui.label(text);
                            });

                            row.col(|ui| {
                                let text = match &data.metrics {
                                    // todo: Cache this calc?
                                    Some(m) => &format!("{}%", m.gc_portion * 100.),
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
        });
    });
}
