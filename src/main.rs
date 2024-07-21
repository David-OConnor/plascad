//! "yeah i would just use AF2 and homology for this"

//! todo: Sort out when to use `bio` types.

// Disables the terminal window. Use this for releases, but disable when debugging.
// #![windows_subsystem = "windows"]

use bincode::{Decode, Encode};
// use bio::{
//     bio_types::sequence::{Sequence, SequenceRead},
//     data_structures::fmindex::FMIndexable,
//     io::fastq::FastqRead,
// };
use eframe::{self, egui, egui::Context};
use primer::PrimerData;

use crate::{gui::PagePrimerCreation, pcr::PolymeraseType, primer::TM_TARGET, util::load};
// use image::GenericImageView;
use crate::{
    gui::{Page, WINDOW_HEIGHT, WINDOW_TITLE, WINDOW_WIDTH},
    pcr::PcrParams,
};

mod gui;
mod primer;
mod solution_helper;
mod util;
// mod snapgene_parse;
mod pcr;
mod toxic_proteins;

// Index 0: 5' end.
type Seq = Vec<Nucleotide>;

/// A DNA nucleotide.
/// todo: RNA A/R
#[derive(Clone, Copy, PartialEq, Debug, Encode, Decode)]
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
    }
}

/// Variables for UI fields, for determining PCR parameters.
#[derive(Encode, Decode)]
struct PcrUi {
    pub primer_tm: f32,
    pub product_len: usize,
    pub polymerase_type: PolymeraseType,
    pub num_cycles: u16,
}

impl Default for PcrUi {
    fn default() -> Self {
        Self {
            primer_tm: TM_TARGET,
            product_len: 1_000,
            polymerase_type: Default::default(),
            num_cycles: 30,
        }
    }
}

#[derive(Default, Encode, Decode)]
struct StateUi {
    // todo: Make separate primer cols and primer data; data in state. primer_cols are pre-formatted
    // todo to save computation.
    page: Page,
    page_primer_creation: PagePrimerCreation,
    seq_insert_input: String,
    seq_vector_input: String,
    seq_amplicon_input: String,
    pcr: PcrUi,
    // pcr_primer: Option<usize>, // primer index, if primer count > 0.
    pcr_primer: usize, // primer index
    primer_selected: Option<usize>,
}

/// Note: use of serde traits here and on various sub-structs are for saving and loading.
#[derive(Default, Encode, Decode)]
struct State {
    ui: StateUi,
    primer_data: Vec<PrimerData>,
    /// Insert and vector are for SLIC and FC.
    seq_insert: Seq,
    seq_vector: Seq,
    /// Amplicon is for basic PCR.
    seq_amplicon: Seq,
    insert_loc: usize,
    /// These limits for choosing the insert location may be defined by the vector's promoter, RBS etc.
    insert_location_5p_limit: usize,
    insert_location_3p_limit: usize,
    pcr: PcrParams,
}

impl State {
    /// Update sequences based on input strings.
    /// todo: This fn may no longer be necessary.
    pub fn _sync_seqs(&mut self) {
        self.seq_insert = util::seq_from_str(&self.ui.seq_insert_input);
        self.seq_vector = util::seq_from_str(&self.ui.seq_vector_input);
        self.seq_amplicon = util::seq_from_str(&self.ui.seq_amplicon_input);
    }

    /// Runs the match serach between primers and sequences. Run this when primers and sequences change.
    pub fn sync_primer_matches(&mut self, primer_i: Option<usize>) {
        let p_list = match primer_i {
            Some(i) => &mut self.primer_data[i..i+1],
            None => &mut self.primer_data,
        };

        for p_data in p_list {
            p_data.matches_amplification_seq = p_data.primer.match_to_seq(&self.seq_amplicon);
            p_data.matches_slic_insert = p_data.primer.match_to_seq(&self.seq_insert);
            p_data.matches_slic_vector = p_data.primer.match_to_seq(&self.seq_vector);
        }
    }

    pub fn sync_pcr(&mut self) {
        self.pcr = PcrParams::new(&self.ui.pcr);
    }

    pub fn sync_metrics(&mut self) {
        for primer in &mut self.primer_data {
            if let Some(metrics) = &mut primer.metrics {
                metrics.update_scores();
            }
        }
    }
}

fn main() {
    // todo: Move to a more robust save/load system later.

    let mut state = load("plasmid_tools.save").unwrap_or_else(|_| State::default());

    // state.sync_seqs();
    state.sync_pcr();
    state.sync_metrics();

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
