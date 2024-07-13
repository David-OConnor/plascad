//! "yeah i would just use AF2 and homology for this"

//! todo: Sort out when to use `bio` types.

mod snapgene_parse;
mod toxic_proteins;

use bio::{
    bio_types::sequence::{Sequence, SequenceRead},
    data_structures::fmindex::FMIndexable,
    io::fastq::FastqRead,
};
use eframe::{self, egui};

use crate::gui::{WINDOW_HEIGHT, WINDOW_TITLE, WINDOW_WIDTH};

mod gui;
mod primer;

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
    seq_expected: Sequence,
    seq_read_assembled: Sequence,
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

struct State {}

fn main() {
    let mut state = State::default();

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([WINDOW_WIDTH, WINDOW_HEIGHT]),
        // icon: load_icon(Path::new("../resources/icon.png")),
        // icon_data,
        follow_system_theme: false,
        ..Default::default()
    };

    eframe::run_native(WINDOW_TITLE, options, Box::new(|_cc| Box::new(state))).unwrap();
}
