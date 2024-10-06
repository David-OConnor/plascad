use bio::{
    alignment::{
        distance::simd::{bounded_levenshtein, hamming, levenshtein},
        pairwise::Aligner,
        Alignment,
    },
    scores::blosum62,
};
use na_seq::{seq_to_letter_bytes, Nucleotide};

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
// pub fn distance(alpha: &[Nucleotide], beta: &[Nucleotide], dist_type: DistanceType) -> u64 {
pub fn distance(alpha: &[Nucleotide], beta: &[Nucleotide]) -> u64 {
    let alpha_ = seq_to_letter_bytes(alpha);
    let beta_ = seq_to_letter_bytes(beta);

    let dist_type = if alpha.len() == beta.len() {
        DistanceType::Hamming
    } else {
        DistanceType::Levenshtein
    };

    match dist_type {
        DistanceType::Hamming => hamming(&alpha_, &beta_),
        DistanceType::Levenshtein => levenshtein(&alpha_, &beta_) as u64,
        DistanceType::BoundedLevenshtein(k) => {
            bounded_levenshtein(&alpha_, &beta_, k).unwrap() as u64
        }
    }
}

// todo: Move some of this to a new `src/alignment` module A/R.
pub fn align_pairwise(seq_0: &[Nucleotide], seq_1: &[Nucleotide]) -> (Alignment, String) {
    // todo: Lots of room to configure this.
    let seq_0_ = seq_to_letter_bytes(seq_0);
    let seq_1_ = seq_to_letter_bytes(seq_1);

    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };

    let mut aligner = Aligner::with_capacity(seq_0.len(), seq_1.len(), -5, -1, &score);

    // todo: Global? Semiglobal? Local?
    let result = aligner.semiglobal(&seq_0_, &seq_1_);

    let text = result.pretty(&seq_0_, &seq_1_, 120);

    (result, text)
}
