//! "yeah i would just use AF2 and homology for this"

//! todo: Sort out when to use `bio` types.

// Disables the terminal window. Use this for releases, but disable when debugging.
// #![windows_subsystem = "windows"]

use std::{fs, fs::File, io::Read, path::Path};

use eframe::egui::IconData;
// use bio::{
//     bio_types::sequence::{Sequence, SequenceRead},
//     data_structures::fmindex::FMIndexable,
//     io::fastq::FastqRead,
// };
use eframe::{self, egui, egui::Context};
use image::{GenericImageView, ImageError};

use crate::gui::{Page, PrimerTableCol, WINDOW_HEIGHT, WINDOW_TITLE, WINDOW_WIDTH};

mod gui;
mod primer;
mod solution_helper;
mod util;
// mod snapgene_parse;
mod toxic_proteins;

type Seq = Vec<Nucleotide>;

/// A DNA nucleotide.
/// todo: RNA A/R
#[derive(Clone, Copy, PartialEq)]
enum Nucleotide {
    A,
    T,
    G,
    C,
}

impl Nucleotide {
    pub fn as_str(&self) -> &str {
        match self {
            Self::A => "a",
            Self::T => "t",
            Self::C => "c",
            Self::G => "g",
        }
    }
}

pub fn make_seq_str(seq: &Seq) -> String {
    let mut result = String::new();

    for nt in seq {
        result.push_str(nt.as_str());
    }

    result
}

pub fn seq_from_str(str: &str) -> Seq {
    let mut result = Vec::new();

    for char in str.to_lowercase().chars() {
        match char {
            'a' => result.push(Nucleotide::A),
            't' => result.push(Nucleotide::T),
            'c' => result.push(Nucleotide::C),
            'g' => result.push(Nucleotide::G),
            _ => (),
        };
    }

    result
}

#[derive(Clone, Copy)]
enum AminoAcid {
    Met,
    Ser,
    Ile,
    Gln,
    His,
    Phe,
    Arg,
    Val,
    Leu,
    Cys,
    Asp,
}

#[derive(Clone, Copy)]
enum Enzyme {
    AsiSi,
    BamHI,
    BclI,
    BlpI,
    BmtI,
}

struct PlasmidData {
    // todo : Which seq type? There are many.
    seq_expected: Seq,
    seq_read_assembled: Seq,
    // seq_reads: Vec<SequenceRead>,
}

fn check_seq_integrity(data: &PlasmidData) {}

fn check_toxic_proteins(data: &PlasmidData) {}

fn check_all(data: &PlasmidData) {
    check_seq_integrity(data);
    check_toxic_proteins(data);
}

impl eframe::App for State {
    /// This is the GUI's event loop.
    fn update(&mut self, ctx: &Context, _frame: &mut eframe::Frame) {
        gui::draw(self, ctx);

        // This repaint is required to keep the event loop running, even with no GUI activity.
        ctx.request_repaint()
    }
}

#[derive(Default)]
struct StateUi {
    // todo: Make separate primer cols and primer data; data in state. primer_cols are pre-formatted
    // todo to save computation.
    primer_cols: Vec<PrimerTableCol>,
    page: Page,
}

#[derive(Default)]
struct State {
    ui: StateUi,
    seq_insert: Seq,
    seq_vector: Seq,
    /// These limits for choosing the insert location may be defined by the vector's promoter, RBS etc.
    insert_location_5p_limit: usize,
    insert_location_3p_limit: usize,
}

fn main() {
    let mut state = State::default();

    let icon_bytes: &[u8] = include_bytes!("resources/icon.png");
    let icon_data = eframe::icon_data::from_png_bytes(icon_bytes);

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([WINDOW_WIDTH, WINDOW_HEIGHT])
            .with_icon(icon_data.unwrap()),
        follow_system_theme: false,
        ..Default::default()
    };

    eframe::run_native(WINDOW_TITLE, options, Box::new(|_cc| Ok(Box::new(state)))).unwrap();
}
