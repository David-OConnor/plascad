use na_seq::{AminoAcid, Nucleotide, Seq, seq_aa_to_u8_lower, seq_to_u8_lower};

// Pairwise alignment scoring. A gap of length k scores `GAP_OPEN + k * GAP_EXTEND`.
// todo: Lots of room to configure this.
const SCORE_MATCH: i32 = 1;
const SCORE_MISMATCH: i32 = -1;
const GAP_OPEN: i32 = -5;
const GAP_EXTEND: i32 = -1;

/// Characters per row of the alignment's text display.
const DISPLAY_WIDTH: usize = 120;

/// Lower than any real score, but with room to add penalties to without overflowing.
const SCORE_NONE: i32 = i32::MIN / 4;

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

/// The result of a pairwise alignment.
#[derive(Clone, Debug)]
pub struct Alignment {
    pub score: i32,
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

/// The number of positions at which two sequences of equal length differ.
fn hamming(alpha: &[u8], beta: &[u8]) -> u64 {
    alpha.iter().zip(beta).filter(|(a, b)| a != b).count() as u64
}

/// The minimum number of substitutions, insertions, and deletions to turn one sequence into the
/// other.
fn levenshtein(alpha: &[u8], beta: &[u8]) -> u32 {
    let mut prev: Vec<u32> = (0..=beta.len() as u32).collect();
    let mut cur = vec![0; beta.len() + 1];

    for (i, a) in alpha.iter().enumerate() {
        cur[0] = i as u32 + 1;

        for (j, b) in beta.iter().enumerate() {
            let substitution = prev[j] + (a != b) as u32;
            cur[j + 1] = substitution.min(prev[j + 1] + 1).min(cur[j] + 1);
        }

        std::mem::swap(&mut prev, &mut cur);
    }

    prev[beta.len()]
}

/// The Levenshtein distance, if it's at most `k`.
fn bounded_levenshtein(alpha: &[u8], beta: &[u8], k: u32) -> Option<u32> {
    // The distance is at least the difference in length.
    if alpha.len().abs_diff(beta.len()) > k as usize {
        return None;
    }

    let dist = levenshtein(alpha, beta);
    (dist <= k).then_some(dist)
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
        DistanceType::Hamming => hamming(alpha, beta),
        DistanceType::Levenshtein => levenshtein(alpha, beta) as u64,
        DistanceType::BoundedLevenshtein(k) => {
            bounded_levenshtein(alpha, beta, k).map_or(u64::MAX, u64::from)
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

/// The alignment state a cell's best score ends in.
#[derive(Clone, Copy, PartialEq)]
#[repr(u8)]
enum AlignState {
    /// Letters from both sequences aligned together.
    Match = 0,
    /// A letter from the first sequence against a gap.
    GapSecond = 1,
    /// A letter from the second sequence against a gap.
    GapFirst = 2,
}

impl AlignState {
    fn from_bits(v: u8) -> Self {
        match v & 0b11 {
            0 => Self::Match,
            1 => Self::GapSecond,
            _ => Self::GapFirst,
        }
    }
}

/// The highest of the three, and which it is. Ties go to the earlier state.
fn best_of(m: i32, gap_second: i32, gap_first: i32) -> (i32, AlignState) {
    let mut result = (m, AlignState::Match);

    if gap_second > result.0 {
        result = (gap_second, AlignState::GapSecond);
    }
    if gap_first > result.0 {
        result = (gap_first, AlignState::GapFirst);
    }

    result
}

/// One column of an alignment: a letter from each sequence, or `None` for a gap.
type Column = (Option<u8>, Option<u8>);

/// Semiglobal alignment with affine gap penalties (Gotoh). All of `x` is aligned, against any part
/// of `y`; the parts of `y` before and after that are free. Returns the score, and the columns,
/// including `y`'s unaligned ends.
///
/// Memory is one byte per pair of positions, for the traceback.
fn align_semiglobal(x: &[u8], y: &[u8]) -> (i32, Vec<Column>) {
    let (n, m) = (x.len(), y.len());
    let width = m + 1;

    // Scores for the previous and current rows of each state. Row 0 is the start: any position in
    // `y`, as its prefix is free.
    let mut prev_m = vec![0; width];
    let mut prev_gs = vec![SCORE_NONE; width];
    let mut prev_gf = vec![SCORE_NONE; width];

    let mut cur_m = vec![SCORE_NONE; width];
    let mut cur_gs = vec![SCORE_NONE; width];
    let mut cur_gf = vec![SCORE_NONE; width];

    // For each cell, the state each state's best score came from; two bits per state.
    let mut traceback = vec![0u8; (n + 1) * width];

    for i in 1..=n {
        // No letters of `y` yet; `x`'s letters so far are against gaps.
        cur_m[0] = SCORE_NONE;
        cur_gf[0] = SCORE_NONE;
        let (score, from) = best_of(
            prev_m[0] + GAP_OPEN + GAP_EXTEND,
            prev_gs[0] + GAP_EXTEND,
            SCORE_NONE,
        );
        cur_gs[0] = score;
        traceback[i * width] = (from as u8) << 2;

        for j in 1..=m {
            let score_pair = if x[i - 1] == y[j - 1] {
                SCORE_MATCH
            } else {
                SCORE_MISMATCH
            };

            let (score, from_m) = best_of(prev_m[j - 1], prev_gs[j - 1], prev_gf[j - 1]);
            cur_m[j] = score + score_pair;

            let (score, from_gs) = best_of(
                prev_m[j] + GAP_OPEN + GAP_EXTEND,
                prev_gs[j] + GAP_EXTEND,
                prev_gf[j] + GAP_OPEN + GAP_EXTEND,
            );
            cur_gs[j] = score;

            let (score, from_gf) = best_of(
                cur_m[j - 1] + GAP_OPEN + GAP_EXTEND,
                cur_gs[j - 1] + GAP_OPEN + GAP_EXTEND,
                cur_gf[j - 1] + GAP_EXTEND,
            );
            cur_gf[j] = score;

            traceback[i * width + j] = from_m as u8 | (from_gs as u8) << 2 | (from_gf as u8) << 4;
        }

        std::mem::swap(&mut prev_m, &mut cur_m);
        std::mem::swap(&mut prev_gs, &mut cur_gs);
        std::mem::swap(&mut prev_gf, &mut cur_gf);
    }

    // `prev_*` now hold the last row. `y`'s suffix is free, so the alignment can end anywhere in
    // it. (Ending on a gap in `x` would only add a penalty for that free suffix.)
    let mut end = (SCORE_NONE, 0, AlignState::Match);
    for j in 0..=m {
        let (score, state) = best_of(prev_m[j], prev_gs[j], SCORE_NONE);
        if score > end.0 {
            end = (score, j, state);
        }
    }

    let (score, j_end, mut state) = end;
    let score = if n == 0 { 0 } else { score };

    let mut aligned = Vec::new();
    let (mut i, mut j) = (n, j_end);

    while i > 0 {
        let tb = traceback[i * width + j];

        match state {
            AlignState::Match => {
                aligned.push((Some(x[i - 1]), Some(y[j - 1])));
                state = AlignState::from_bits(tb);
                i -= 1;
                j -= 1;
            }
            AlignState::GapSecond => {
                aligned.push((Some(x[i - 1]), None));
                state = AlignState::from_bits(tb >> 2);
                i -= 1;
            }
            AlignState::GapFirst => {
                aligned.push((None, Some(y[j - 1])));
                state = AlignState::from_bits(tb >> 4);
                j -= 1;
            }
        }
    }

    let mut columns: Vec<Column> = y[..j].iter().map(|&b| (None, Some(b))).collect();
    columns.extend(aligned.into_iter().rev());
    columns.extend(y[j_end..].iter().map(|&b| (None, Some(b))));

    (score, columns)
}

/// Three rows per block: the first sequence, markers, and the second sequence. `|` marks a match,
/// `.` a mismatch, and `-` a gap.
fn alignment_text(columns: &[Column]) -> String {
    let mut result = String::new();

    for block in columns.chunks(DISPLAY_WIDTH) {
        let row = |pick: fn(&Column) -> Option<u8>| -> String {
            block
                .iter()
                .map(|c| pick(c).map_or('-', char::from))
                .collect()
        };

        let markers: String = block
            .iter()
            .map(|c| match c {
                (Some(a), Some(b)) if a == b => '|',
                (Some(_), Some(_)) => '.',
                _ => ' ',
            })
            .collect();

        if !result.is_empty() {
            result.push('\n');
        }
        result.push_str(&format!("{}\n{markers}\n{}\n", row(|c| c.0), row(|c| c.1)));
    }

    result
}

fn align_pairwise(seq_0: &[u8], seq_1: &[u8]) -> (Alignment, String) {
    // todo: Global? Semiglobal? Local?
    let (score, columns) = align_semiglobal(seq_0, seq_1);

    (Alignment { score }, alignment_text(&columns))
}

pub fn align_pairwise_nt(seq_0: &[Nucleotide], seq_1: &[Nucleotide]) -> (Alignment, String) {
    let seq_0_ = seq_to_u8_lower(seq_0);
    let seq_1_ = seq_to_u8_lower(seq_1);

    align_pairwise(&seq_0_, &seq_1_)
}

pub fn align_pairwise_aa(seq_0: &[AminoAcid], seq_1: &[AminoAcid]) -> (Alignment, String) {
    let seq_0_ = seq_aa_to_u8_lower(seq_0);
    let seq_1_ = seq_aa_to_u8_lower(seq_1);

    align_pairwise(&seq_0_, &seq_1_)
}
