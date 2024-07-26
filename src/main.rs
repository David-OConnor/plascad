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

use crate::{
    gui::{Page, PageSeq, WINDOW_HEIGHT, WINDOW_TITLE, WINDOW_WIDTH},
    pcr::{PcrParams, PolymeraseType},
    primer::{PrimerDirection, TM_TARGET},
    restriction_enzyme::{load_re_library, ReMatch, RestrictionEnzyme},
    util::{load, seq_complement, seq_from_str},
};

mod gui;
mod melting_temp_calcs;
mod pcr;
mod primer;
mod primer_metrics;
mod restriction_enzyme;
mod solution_helper;
mod toxic_proteins;
mod util;
// mod snapgene_parse;

// Index 0: 5' end.
type Seq = Vec<Nucleotide>;

/// A DNA nucleotide.
/// todo: RNA A/R
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Encode, Decode)]
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
    /// index from primer data. For storing dropdown state.
    pub primer_selected: usize,
}

impl Default for PcrUi {
    fn default() -> Self {
        Self {
            primer_tm: TM_TARGET,
            product_len: 1_000,
            polymerase_type: Default::default(),
            num_cycles: 30,
            primer_selected: 0,
        }
    }
}

#[derive(Encode, Decode)]
/// Concentrations of common ions in the oglio solution. Affects melting temperature (TM).
/// All values are in milliMolar.
struct IonConcentrations {
    /// Na+ or K+
    pub monovalent: f32,
    /// Mg2+
    pub divalent: f32,
    pub dntp: f32,
    /// Primer concentration, in nM.
    pub primer: f32,
}

impl Default for IonConcentrations {
    fn default() -> Self {
        // todo: Adjust A/R
        Self {
            monovalent: 50.,
            divalent: 1.5,
            dntp: 0.2,
            primer: 25.,
        }
    }
}

/// Values defined here generally aren't worth saving to file etc.
#[derive(Default, Encode, Decode)]
struct StateUi {
    // todo: Make separate primer cols and primer data; data in state. primer_cols are pre-formatted
    // todo to save computation.
    page: Page,
    // page_primer: PagePrimer,
    page_seq: PageSeq,
    seq_insert_input: String,
    seq_vector_input: String,
    seq_input: String,
    pcr: PcrUi,
    // pcr_primer: Option<usize>, // primer index, if primer count > 0.
    pcr_primer: usize, // primer index
    primer_selected: Option<usize>,
    // todo: suitstruct for show/hide A/R
    /// Hide the primer addition/creation/QC panel.
    hide_primer_table: bool,
    /// Show or hide restriction enzymes from the sequence view.
    show_res: bool,
    /// Show and hide primers on
    show_primers: bool,
}

/// Note: use of serde traits here and on various sub-structs are for saving and loading.
#[derive(Default, Encode, Decode)]
struct State {
    ui: StateUi, // Does not need to be saved
    primer_data: Vec<PrimerData>,
    /// Insert and vector are for SLIC and FC.
    // seq_insert: Seq,
    // seq_vector: Seq,
    seq: Seq,
    // /// Amplicon is for basic PCR.
    // seq_amplicon: Seq,
    insert_loc: usize,
    /// These limits for choosing the insert location may be defined by the vector's promoter, RBS etc.
    insert_location_5p_limit: usize,
    insert_location_3p_limit: usize,
    ion_concentrations: IonConcentrations,
    pcr: PcrParams,
    restriction_enzyme_lib: Vec<RestrictionEnzyme>, // Does not need to be saved
    restriction_enzyme_sites: Vec<ReMatch>,
}

impl State {
    /// Runs the match serach between primers and sequences. Run this when primers and sequences change.
    pub fn sync_primer_matches(&mut self, primer_i: Option<usize>) {
        let p_list = match primer_i {
            Some(i) => &mut self.primer_data[i..i + 1],
            None => &mut self.primer_data,
        };

        for p_data in p_list {
            p_data.matches_seq = p_data.primer.match_to_seq(&self.seq);
            // p_data.matches_insert = p_data.primer.match_to_seq(&self.seq_insert);
            // p_data.matches_vector = p_data.primer.match_to_seq(&self.seq_vector);
            // p_data.matches_vector_with_insert =
            //     p_data.primer.match_to_seq(&self.seq);
        }
    }

    pub fn sync_pcr(&mut self) {
        self.pcr = PcrParams::new(&self.ui.pcr);
    }

    /// Identify restriction enzyme sites in the sequence
    pub fn sync_re_sites(&mut self) {
        self.restriction_enzyme_sites = Vec::new();

        // let seq_comp = seq_complement(seq);

        for (lib_index, re) in self.restriction_enzyme_lib.iter().enumerate() {
            // todo: Use bio lib?
            for i in 0..self.seq.len() {
                if i + re.seq.len() + 1 >= self.seq.len() {
                    continue;
                }

                if re.seq == self.seq[i..i + re.seq.len()] {
                    self.restriction_enzyme_sites.push(ReMatch {
                        lib_index,
                        seq_index: i,
                        direction: PrimerDirection::Forward,
                    });
                }
                // todo: Simpler way?
                // todo: Evaluate if this is what you want.
                // if re.seq == seq_comp[i..i + re.seq.len()] {
                //     self.restriction_enzyme_sites.push(ReMatch {
                //         lib_index: i_re,
                //         seq_index: i,
                //         direction: PrimerDirection::Reverse,
                //     });
                // }
            }
        }

        // This sorting aids in our up/down label alternation in the display.
        self.restriction_enzyme_sites
            .sort_by(|a, b| a.seq_index.cmp(&b.seq_index));
    }

    pub fn sync_metrics(&mut self) {
        for primer in &mut self.primer_data {
            if let Some(metrics) = &mut primer.metrics {
                metrics.update_scores();
            }
        }
    }

    /// Update the combined SLIC vector + insert sequence.
    pub fn sync_cloning_product(&mut self) {
        let seq_vector = seq_from_str(&self.ui.seq_vector_input);
        let seq_insert = seq_from_str(&self.ui.seq_insert_input);

        self.seq = seq_vector.clone(); // Clone the original vector
        self.seq
            .splice(self.insert_loc..self.insert_loc, seq_insert.iter().cloned());

        self.sync_primer_matches(None);

        //
        // self.seq_vector_with_insert = self.seq_vector[..self.insert_loc].clone();
        //
        // for nt in self.seq_insert {
        //     self.seq_vector_with_insert.push(nt);
        // }
        //
        // for nt in self.seq_vector[self.insert_loc..] {
        //     self.seq_vector_with_insert.push(nt);
        // }
    }
}

fn main() {
    // todo: Move to a more robust save/load system later.

    let mut state = load("plasmid_tools.save").unwrap_or_else(|_| State::default());

    state.restriction_enzyme_lib = load_re_library();

    state.sync_pcr();
    state.sync_metrics();
    state.sync_re_sites();
    // state.sync_cloning_product();

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
