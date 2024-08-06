//! This module contains code related to primer (oglionucleotide) design and QC.

use std::ops::Range;

use bincode::{Decode, Encode};

use crate::{
    gui::primer_qc::DEFAULT_TRIM_AMT,
    primer_metrics::PrimerMetrics,
    sequence::{
        seq_complement, seq_from_str, seq_to_str, Nucleotide,
        Nucleotide::{C, G},
        Seq,
    },
    IonConcentrations, State,
};

// If a primer length is below this, many calculations will be disabled for it.
pub const MIN_PRIMER_LEN: usize = 10;
pub const TM_TARGET: f32 = 59.; // Also used as a default for PCR GUI.

// todo: Sort out your types.

/// These are also relevant for FastCloning.
pub struct SlicPrimers {
    pub vector_fwd: Primer,
    pub vector_rev: Primer,
    pub insert_fwd: Primer,
    pub insert_rev: Primer,
}

/// These are also relevant for FastCloning.
pub struct AmplificationPrimers {
    pub fwd: Primer,
    pub rev: Primer,
}

#[derive(Clone, Copy, Debug, PartialEq, Encode, Decode)]
pub enum PrimerDirection {
    Forward,
    Reverse,
}

#[derive(Default, Clone, Encode, Decode)]
pub struct Primer {
    pub sequence: Seq,
    pub description: String,
}

impl Primer {
    /// Match this primer to a sequence. Check both directions.
    /// Returns direction, and start and end indexes of the sequence. If direction is reversed,
    /// the indexs matches to the reversed index.
    pub fn match_to_seq(&self, seq: &[Nucleotide]) -> Vec<(PrimerDirection, Range<usize>)> {
        let mut result = Vec::new();

        // This check prevents spurious small-sequence matches, which may be numerous otherwise.
        if self.sequence.len() < MIN_PRIMER_LEN {
            return result;
        }

        // todo: Partial matches as well.
        let seq_len = seq.len();
        let primer_len = self.sequence.len();
        let complement = seq_complement(seq);

        for seq_start in 0..seq_len {
            // Note: This approach handles sequence wraps, eg [circular] plasmids.
            let seq_iter = seq.iter().cycle().skip(seq_start).take(primer_len);
            if self.sequence.iter().eq(seq_iter) {
                let seq_end = (seq_start + primer_len) % seq_len;
                result.push((PrimerDirection::Forward, seq_start..seq_end));
            }

            // let end_i = (seq_start_i + self.sequence.len()) % seq_len;
            // if self.sequence == seq[seq_start_i..end_i] {
            //     result.push((PrimerDirection::Forward, seq_start_i))
            // }
        }

        for seq_start in 0..seq_len {
            let seq_iter = complement.iter().cycle().skip(seq_start).take(primer_len);
            if self.sequence.iter().eq(seq_iter) {
                let seq_end = (seq_start + primer_len) % seq_len;
                result.push((PrimerDirection::Reverse, seq_start..seq_end));
            }
        }

        result
    }
}

pub fn design_slic_fc_primers(
    seq_vector: &Seq,
    seq_insert: &Seq,
    insert_loc: usize,
) -> Option<SlicPrimers> {
    // These lenghts should be long enough for reasonablely high-length primers, should that be
    // required for optimal characteristics.
    const UNTRIMMED_LEN_INSERT: usize = 30;
    const UNTRIMMED_LEN_VECTOR: usize = 32;

    let seq_len_vector = seq_vector.len();
    let seq_len_insert = seq_insert.len();

    if insert_loc > seq_len_vector {
        eprintln!("Invalid insert loc. Loc: {insert_loc}, vector: {seq_len_vector}");
        // todo: Return an error, and show that in the UI.
        return None;
    }

    let vector_reversed = seq_complement(seq_vector);
    let insert_reversed = seq_complement(seq_insert);

    // todo: You should wrap the vector if the insert loc is near the edges. Clamping
    // todo as we do will produce incorrect behavior.

    let insert_loc_reversed = seq_len_vector - insert_loc;

    let (seq_vector_fwd, seq_vector_rev) = {
        let mut end = insert_loc + UNTRIMMED_LEN_VECTOR;
        end = end.clamp(0, seq_len_vector);

        let mut end_reversed = insert_loc_reversed + UNTRIMMED_LEN_VECTOR;
        end_reversed = end_reversed.clamp(0, seq_len_vector);

        (
            seq_vector[insert_loc..end].to_owned(),
            vector_reversed[insert_loc_reversed..end_reversed].to_owned(),
        )
    };

    let (seq_insert_fwd, seq_insert_rev) = {
        let mut insert_end = UNTRIMMED_LEN_INSERT;
        insert_end = insert_end.clamp(0, seq_len_insert);

        let mut insert_end_reversed = UNTRIMMED_LEN_INSERT;
        insert_end_reversed = insert_end_reversed.clamp(0, seq_len_insert);

        // We will combine these with vector seqs for the final insert primers.
        let insert_only_seq_fwd = &seq_insert[..insert_end];
        let insert_only_seq_rev = &insert_reversed[..insert_end_reversed];

        let mut fwd = seq_complement(&seq_vector_rev);
        fwd.extend(insert_only_seq_fwd);

        let mut rev = seq_complement(&seq_vector_fwd);
        rev.extend(insert_only_seq_rev);

        (fwd, rev)
    };

    // We tune downstream.
    Some(SlicPrimers {
        vector_fwd: Primer {
            sequence: seq_vector_fwd,
            description: "Vector fwd".to_owned(),
        },
        vector_rev: Primer {
            sequence: seq_vector_rev,
            description: "Vector rev".to_owned(),
        },
        insert_fwd: Primer {
            sequence: seq_insert_fwd,
            description: "Insert fwd".to_owned(),
        },
        insert_rev: Primer {
            sequence: seq_insert_rev,
            description: "Insert rev".to_owned(),
        },
    })
}

// todo: Use this A/R, called from the UI page.
pub fn design_amplification_primers(seq: &[Nucleotide]) -> Option<AmplificationPrimers> {
    // These lenghts should be long enough for reasonablely high-length primers, should that be
    // required for optimal characteristics.
    const UNTRIMMED_LEN: usize = 32;

    let seq_len = seq.len();
    let reversed = seq_complement(seq);

    let (seq_fwd, seq_rev) = {
        let mut end = UNTRIMMED_LEN;
        end = end.clamp(0, seq_len);

        let mut end_reversed = UNTRIMMED_LEN;
        end_reversed = end_reversed.clamp(0, seq_len);

        (seq[..end].to_owned(), reversed[..end_reversed].to_owned())
    };

    // We tune downstream.
    Some(AmplificationPrimers {
        fwd: Primer {
            sequence: seq_fwd,
            description: "Amplification fwd".to_owned(),
        },
        rev: Primer {
            sequence: seq_rev,
            description: "Amplification rev".to_owned(),
        },
    })
}

#[derive(Default, Clone, Encode, Decode)]
pub struct PrimerData {
    /// This primer excludes nts past the tuning offsets.
    pub primer: Primer,
    /// Editing is handled using this string; we convert the string to our nucleotide sequence as needed.
    /// This includes nts past the tuning offsets.
    pub sequence_input: String,
    pub metrics: Option<PrimerMetrics>,
    // tunable_end: TunableEnd,
    /// These fields control if a given primer end is fixed (Eg marking the start of an insert,
    /// marking the insert point in a vector etc) or if we can tune its length to optimize the primer.
    pub tunable_5p: TuneSetting,
    pub tunable_3p: TuneSetting,
    /// These seq_removed fields are redundant with primer and tune settings. We use them to cache
    /// the actual sequences that are removed for display purposes.
    // seq_removed_5p: Seq,
    // seq_removed_3p: Seq
    pub seq_removed_5p: String,
    pub seq_removed_3p: String,
    /// todo: Which direction is the range, if the direction is reverse?
    pub matches_seq: Vec<(PrimerDirection, Range<usize>)>,
    // pub matches_vector: Vec<(PrimerDirection, Range<usize>)>, // todo: Currently unused.
    // pub matches_insert: Vec<(PrimerDirection, Range<usize>)>, // todo: Currently unused.
    // pub matches_vector_with_insert: Vec<(PrimerDirection, Range<usize>)>,
}

impl PrimerData {
    pub fn new(primer: Primer) -> Self {
        let mut result = Self::default();
        result.primer = primer;
        result.sequence_input = seq_to_str(&result.primer.sequence);

        result
    }

    /// Perform calculations on primer quality and related data. Run this when the sequence changes,
    /// the tuning values change etc.
    pub fn run_calcs(&mut self, ion_concentrations: &IonConcentrations) {
        let full_len = self.sequence_input.len();
        let mut start = 0;
        let mut end = full_len;

        if let TuneSetting::Enabled(i) = self.tunable_5p {
            start = i;
        }

        if let TuneSetting::Enabled(i) = self.tunable_3p {
            end = if i > full_len {
                // Prevents an overrun.
                0
            } else {
                full_len - i
            };
        }

        if start > end || start + 1 > self.sequence_input.len() {
            start = 0;
            end = full_len
        }

        self.primer.sequence = seq_from_str(&self.sequence_input[start..end]);
        self.metrics = self.primer.calc_metrics(ion_concentrations);

        self.sequence_input[..start].clone_into(&mut self.seq_removed_5p);
        self.sequence_input[end..].clone_into(&mut self.seq_removed_3p);
    }

    /// Automatically select primer length based on quality score.
    pub fn tune(&mut self, ion: &IonConcentrations) {
        if self.tunable_3p != TuneSetting::Disabled && self.tunable_5p == TuneSetting::Disabled {
            self.tune_single_end(ion);
        }
        if self.tunable_3p == TuneSetting::Disabled && self.tunable_5p != TuneSetting::Disabled {
            self.tune_single_end(ion);
        }
        if self.tunable_3p != TuneSetting::Disabled && self.tunable_5p != TuneSetting::Disabled {
            self.tune_both_ends(ion);
        }
    }

    fn tune_single_end(&mut self, ion: &IonConcentrations) {
        let primer_len = self.primer.sequence.len();
        if primer_len <= MIN_PRIMER_LEN {
            return;
        }

        let mut best_val = 0;
        let mut best_score = 0.;

        let num_vals = primer_len - MIN_PRIMER_LEN;

        for val in 0..num_vals {
            // When this function is called, exactly one of these ends should be enabled.
            if let TuneSetting::Enabled(tune_val) = &mut self.tunable_3p {
                *tune_val = val;
            }
            if let TuneSetting::Enabled(tune_val) = &mut self.tunable_5p {
                *tune_val = val;
            }
            self.run_calcs(ion);

            if let Some(metrics) = &self.metrics {
                if metrics.quality_score > best_score {
                    best_val = val;
                    best_score = metrics.quality_score;
                }
            }
        }

        if let TuneSetting::Enabled(tune_val) = &mut self.tunable_3p {
            // This is getting a bit verbose.
            *tune_val = best_val;
        }

        self.run_calcs(ion);
    }

    fn tune_both_ends(&mut self, ion: &IonConcentrations) {
        let primer_len = self.primer.sequence.len();
        if primer_len <= MIN_PRIMER_LEN {
            return;
        }

        let mut best_val_5p = 0;
        let mut best_val_3p = 0;
        let mut best_score = 0.;

        // We will iterate over all possible tunings within the bounds of the limits of both ends,
        // and minimum primer length. We do so by varying the start point (We set this, starting at
        // the 5' end and expanding towards 3', but either direction will do), and length, up
        // to the other end.

        for val_5p in 0..primer_len - MIN_PRIMER_LEN {
            for len in MIN_PRIMER_LEN..primer_len - val_5p {
                let val_3p = primer_len - (val_5p + len);

                if let TuneSetting::Enabled(tune_val) = &mut self.tunable_3p {
                    *tune_val = val_3p;
                }
                if let TuneSetting::Enabled(tune_val) = &mut self.tunable_5p {
                    *tune_val = val_5p;
                }
                self.run_calcs(ion);

                if let Some(metrics) = &self.metrics {
                    if metrics.quality_score > best_score {
                        best_val_5p = val_5p;
                        best_val_3p = val_3p;
                        best_score = metrics.quality_score;
                    }
                }
            }
        }

        if let TuneSetting::Enabled(tune_val) = &mut self.tunable_5p {
            // This is getting a bit verbose.
            *tune_val = best_val_5p;
        }

        if let TuneSetting::Enabled(tune_val) = &mut self.tunable_3p {
            // This is getting a bit verbose.
            *tune_val = best_val_3p;
        }
        self.run_calcs(ion);
    }
}

impl Default for TuneSetting {
    fn default() -> Self {
        Self::Disabled
    }
}

#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum TuneSetting {
    Disabled,
    /// Inner: Offset index; this marks the distance from the respective ends that the sequence is attentuated to.
    Enabled(usize),
}

impl TuneSetting {
    pub fn toggle(&mut self) {
        *self = match self {
            Self::Disabled => Self::Enabled(0),
            _ => Self::Disabled,
        }
    }
}

/// Calculate GC portion, on a scale of 0 to 1.
/// This is a standalone fn, as it's used outside of Primer methods.
pub fn calc_gc(seq: &[Nucleotide]) -> f32 {
    let mut num_gc = 0;
    for nt in seq {
        if *nt == C || *nt == G {
            num_gc += 1;
        }
    }

    num_gc as f32 / seq.len() as f32
}

/// We run this to generate cloning primers when clicking the button
pub fn make_cloning_primers(state: &mut State) {
    let seq_vector = seq_from_str(&state.ui.seq_vector_input);
    let seq_insert = seq_from_str(&state.ui.seq_insert_input);

    if let Some(primers) = design_slic_fc_primers(&seq_vector, &seq_insert, state.insert_loc) {
        let sequence_input = seq_to_str(&primers.insert_fwd.sequence);

        let mut insert_fwd = PrimerData {
            primer: primers.insert_fwd,
            sequence_input,
            // Both ends are  tunable, since this glues the insert to the vector
            tunable_5p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
            tunable_3p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
            ..Default::default()
        };

        let sequence_input = seq_to_str(&primers.insert_rev.sequence);
        let mut insert_rev = PrimerData {
            primer: primers.insert_rev,
            sequence_input,
            // Both ends are tunable, since this glues the insert to the vector
            tunable_5p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
            tunable_3p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
            ..Default::default()
        };

        let sequence_input = seq_to_str(&primers.vector_fwd.sequence);
        let mut vector_fwd = PrimerData {
            primer: primers.vector_fwd,
            sequence_input,
            // 5' is non-tunable: This is the insert location.
            tunable_5p: TuneSetting::Disabled,
            tunable_3p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
            ..Default::default()
        };

        let sequence_input = seq_to_str(&primers.vector_rev.sequence);
        let mut vector_rev = PrimerData {
            primer: primers.vector_rev,
            sequence_input,
            tunable_5p: TuneSetting::Disabled,
            // 3' is non-tunable: This is the insert location.
            tunable_3p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
            ..Default::default()
        };

        insert_fwd.tune(&state.ion_concentrations);
        insert_rev.tune(&state.ion_concentrations);
        vector_fwd.tune(&state.ion_concentrations);
        vector_rev.tune(&state.ion_concentrations);

        state
            .primer_data
            .extend([insert_fwd, insert_rev, vector_fwd, vector_rev]);

        state.sync_primer_matches(None); // note: Not requried to run on all primers.
    }
}

pub fn make_amplification_primers(state: &mut State) {
    if let Some(primers) = design_amplification_primers(&state.seq) {
        let sequence_input = seq_to_str(&primers.fwd.sequence);

        let mut primer_fwd = PrimerData {
            primer: primers.fwd,
            sequence_input,
            tunable_5p: TuneSetting::Disabled,
            tunable_3p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
            ..Default::default()
        };

        let sequence_input = seq_to_str(&primers.rev.sequence);
        let mut primer_rev = PrimerData {
            primer: primers.rev,
            sequence_input,
            tunable_5p: TuneSetting::Disabled,
            tunable_3p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
            ..Default::default()
        };

        primer_fwd.tune(&state.ion_concentrations);
        primer_rev.tune(&state.ion_concentrations);

        state.primer_data.extend([primer_fwd, primer_rev]);

        state.sync_primer_matches(None); // note: Not requried to run on all primers.
    }
}
