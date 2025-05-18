use bio::alignment::{
    Alignment,
    distance::simd::{bounded_levenshtein, hamming, levenshtein},
    pairwise::Aligner,
};
use na_seq::{AminoAcid, Nucleotide, Seq, seq_aa_to_u8_lower, seq_to_u8_lower};

#[derive(Clone, Copy, PartialEq)]
pub enum AlignmentMode {
    Dna,
    AminoAcid,
}

impl Default for AlignmentMode {
    fn default() -> Self {
        Self::Dna
    }
}

#[derive(Default)]
pub struct AlignmentState {
    pub seq_a: Seq,
    pub seq_b: Seq,
    pub seq_aa_a: Vec<AminoAcid>,
    pub seq_aa_b: Vec<AminoAcid>,
    // todo: Perhaps these inputs are better fit for State_ui.
    pub seq_a_input: String,
    pub seq_b_input: String,
    pub alignment_result: Option<Alignment>,
    pub dist_result: Option<u64>,
    pub text_display: String, // Ie `AlignmentResult::pretty`.
    pub mode: AlignmentMode,
}

#[derive(Clone, Copy)]
/// Use Hamming Distance if:
///     Sequences are of the same length.
///     You're only looking for substitutions and the sequences are already aligned.
///     Fast, low-complexity calculations are a priority.
/// Use Levenshtein Distance if:
///     Sequences vary in length.
///     Insertions, deletions, and substitutions are all meaningful for your comparison.
///     You need a precise measure of distance for downstream processing.
/// Use Bounded Levenshtein Distance if:
///     You have a known tolerance or threshold for "closeness."
///     You want to filter out sequences that are too different without computing the full distance.
///     Speed and scalability are critical considerations.
pub enum DistanceType {
    Hamming,
    Levenshtein,
    /// Inner: u32
    BoundedLevenshtein(u32),
}

// todo: Which dist algo? bounded_levenshtein, hamming, levenshtein?
/// Accepts UTF-8 byte representation of letters.
fn distance(alpha: &[u8], beta: &[u8]) -> u64 {
    let dist_type = if alpha.len() == beta.len() {
        DistanceType::Hamming
    } else {
        DistanceType::Levenshtein
    };

    match dist_type {
        DistanceType::Hamming => hamming(&alpha, &beta),
        DistanceType::Levenshtein => levenshtein(&alpha, &beta) as u64,
        DistanceType::BoundedLevenshtein(k) => {
            bounded_levenshtein(&alpha, &beta, k).unwrap() as u64
        }
    }
}

pub fn distance_nt(alpha: &[Nucleotide], beta: &[Nucleotide]) -> u64 {
    let alpha_ = seq_to_u8_lower(alpha);
    let beta_ = seq_to_u8_lower(beta);
    distance(&alpha_, &beta_)
}

pub fn distance_aa(alpha: &[AminoAcid], beta: &[AminoAcid]) -> u64 {
    let alpha_ = seq_aa_to_u8_lower(alpha);
    let beta_ = seq_aa_to_u8_lower(beta);
    distance(&alpha_, &beta_)
}

fn align_pairwise(seq_0: &[u8], seq_1: &[u8]) -> (Alignment, String) {
    // todo: Lots of room to configure this.
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };

    let mut aligner = Aligner::with_capacity(seq_0.len(), seq_1.len(), -5, -1, &score);

    // todo: Global? Semiglobal? Local?
    let result = aligner.semiglobal(seq_0, seq_1);

    let text = result.pretty(seq_0, seq_1, 120);

    (result, text)
}

pub fn align_pairwise_nt(seq_0: &[Nucleotide], seq_1: &[Nucleotide]) -> (Alignment, String) {
    // todo: Lots of room to configure this.
    let seq_0_ = seq_to_u8_lower(seq_0);
    let seq_1_ = seq_to_u8_lower(seq_1);

    align_pairwise(&seq_0_, &seq_1_)
}

pub fn align_pairwise_aa(seq_0: &[AminoAcid], seq_1: &[AminoAcid]) -> (Alignment, String) {
    // todo: Lots of room to configure this.
    let seq_0_ = seq_aa_to_u8_lower(seq_0);
    let seq_1_ = seq_aa_to_u8_lower(seq_1);

    align_pairwise(&seq_0_, &seq_1_)
}
