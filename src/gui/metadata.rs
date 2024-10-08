//! Contains references, comments, etc about the plasmid.

use chrono::NaiveDate;
use eframe::egui::{Color32, RichText, ScrollArea, TextEdit, Ui};

use crate::{
    gui::{COL_SPACING, ROW_SPACING},
    misc_types::Metadata,
};

const WIDTH_RATIO: f32 = 0.6;
const ROW_HEIGHT: usize = 1;

const HEADING_COLOR: Color32 = Color32::from_rgb(40, 180, 255);

/// A convenience function to create a text edit for Option<String>
fn option_edit(
    val: &mut Option<String>,
    label: &str,
    multi: bool,
    width_: Option<f32>,
    ui: &mut Ui,
) {
    ui.horizontal(|ui| {
        // ui.allocate_exact_size(Vec2::new(LABEL_WIDTH, 0.0), egui::Sense::hover()); // Reserve space
        ui.label(label);

        let width = match width_ {
            Some(w) => w,
            None => ui.available_width() * WIDTH_RATIO,
        };

        // todo: Way without cloning?
        let mut v = val.clone().unwrap_or_default();
        // Don't use these margins if there is a narrow window.
        let response = if multi {
            ui.add(
                TextEdit::multiline(&mut v)
                    .desired_width(width)
                    .desired_rows(ROW_HEIGHT),
            )
        } else {
            ui.add(TextEdit::singleline(&mut v).desired_width(width))
        };

        if response.changed() {
            *val = if !v.is_empty() {
                Some(v.to_owned())
            } else {
                None
            };
        }
    });

    ui.add_space(ROW_SPACING / 2.);
}

pub fn metadata_page(data: &mut Metadata, ui: &mut Ui) {
    // todo: YOu need neat grid alignment. How can we make the labels take up constant space?

    // todo: Examine which fields should be single vs multiline, and the order.
    ScrollArea::vertical().show(ui, |ui| {
        ui.heading(RichText::new("General:").color(HEADING_COLOR));
        ui.add_space(ROW_SPACING / 2.);

        ui.horizontal(|ui| {
            ui.label("Plasmid name:");
            ui.text_edit_singleline(&mut data.plasmid_name);
        });
        ui.add_space(ROW_SPACING);

        option_edit(&mut data.definition, "Definition:", false, None, ui);

        ui.horizontal(|ui| {
            ui.label("Division:");
            ui.add(TextEdit::singleline(&mut data.division).desired_width(60.));
            ui.add_space(COL_SPACING);

            option_edit(
                &mut data.molecule_type,
                "Molecule type:",
                false,
                Some(80.),
                ui,
            );

            if let Some(date) = data.date {
                ui.add_space(COL_SPACING);
                // todo: Allow the date to be editable, and not just when present.
                let date = NaiveDate::from_ymd(date.0, date.1.into(), date.2.into());
                ui.label(format!("Date: {}", date));
            }
        });
        ui.add_space(ROW_SPACING / 2.);

        option_edit(&mut data.accession, "Accession:", true, None, ui);
        option_edit(&mut data.version, "Version:", true, None, ui);
        option_edit(&mut data.keywords, "Keywords:", true, None, ui);
        option_edit(&mut data.source, "Source:", true, None, ui);
        option_edit(&mut data.organism, "Organism:", true, None, ui);

        ui.add_space(ROW_SPACING);

        //  pub locus: String,
        //     pub definition: Option<String>,
        //     pub accession: Option<String>,
        //     pub version: Option<String>,
        //     // pub keywords: Vec<String>,
        //     pub keywords: Option<String>, // todo vec?
        //     pub source: Option<String>,
        //     pub organism: Option<String>,

        ui.heading(RichText::new("References:").color(HEADING_COLOR));
        ui.add_space(ROW_SPACING / 2.);

        for ref_ in &mut data.references {
            ui.horizontal(|ui| {
                ui.label("Title:");
                ui.add(
                    TextEdit::multiline(&mut ref_.title)
                        .desired_width(ui.available_width() * WIDTH_RATIO)
                        .desired_rows(ROW_HEIGHT),
                );
            });
            ui.add_space(ROW_SPACING / 2.);

            ui.horizontal(|ui| {
                ui.label("Description:");
                ui.add(
                    TextEdit::multiline(&mut ref_.description)
                        .desired_width(ui.available_width() * WIDTH_RATIO)
                        .desired_rows(ROW_HEIGHT),
                );
            });
            ui.add_space(ROW_SPACING / 2.);

            option_edit(&mut ref_.authors, "Authors:", true, None, ui);
            option_edit(&mut ref_.consortium, "Consortium:", true, None, ui);

            option_edit(&mut ref_.journal, "Journal:", true, None, ui);
            option_edit(&mut ref_.pubmed, "Pubmed:", true, None, ui);
            option_edit(&mut ref_.remark, "Remarks:", true, None, ui);

            ui.add_space(ROW_SPACING);
        }

        // egui::Shape::hline(2, 2., Stroke::new(2., Color32::WHITE));

        ui.heading(RichText::new("Comments:").color(HEADING_COLOR));
        ui.add_space(ROW_SPACING);
        if ui.button("➕ Add").clicked() {
            data.comments.push(String::new());
        }

        for comment in &mut data.comments {
            let response = ui.add(
                TextEdit::multiline(comment)
                    .desired_width(ui.available_width() * WIDTH_RATIO)
                    .desired_rows(ROW_HEIGHT),
            );
            // if response.changed() {
            // }

            ui.add_space(ROW_SPACING / 2.);
        }
    });
}
