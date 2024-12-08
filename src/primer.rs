//! This module contains code related to primer (oglionucleotide) design and QC.

use bincode::{Decode, Encode};
use eframe::egui::Color32;
use na_seq::{seq_complement, seq_from_str, seq_to_str_lower, seq_weight, Nucleotide, Seq};

use crate::{
    gui::{primer_table::DEFAULT_TRIM_AMT, PRIMER_FWD_COLOR, PRIMER_REV_COLOR},
    primer_metrics::PrimerMetrics,
    util::{match_subseq, RangeIncl},
    State,
};

// If a primer length is below this, many calculations will be disabled for it.
pub const MIN_PRIMER_LEN: usize = 10;
pub const TM_TARGET: f32 = 59.; // Also used as a default for PCR GUI.

// These lenghts should be long enough for reasonablely high-length primers, should that be
// required for optimal characteristics.
const UNTRIMMED_LEN_INSERT: usize = 30;
const UNTRIMMED_LEN_VECTOR: usize = 32;

// todo: Sort out your types.

#[derive(Clone, Debug, Encode, Decode)]
pub struct PrimerMatch {
    pub direction: PrimerDirection,
    /// This range is in the forward direction for both primers.
    pub range: RangeIncl,
}

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

impl PrimerDirection {
    pub fn color(&self) -> Color32 {
        match self {
            Self::Forward => PRIMER_FWD_COLOR,
            Self::Reverse => PRIMER_REV_COLOR,
        }
    }
}

#[derive(Default, Clone, Encode, Decode)]
pub struct Primer {
    pub sequence: Seq,
    pub name: String,
    pub description: Option<String>, // todo: Display this.
    /// Data that is dynamically regenerated, and generally not important for saving and loading to files.
    pub volatile: PrimerData,
}

impl Primer {
    /// Match this primer to a sequence. Check both directions.
    /// Returns direction, and start and end indexes of the sequence. If direction is reversed,
    /// the indexs matches to the reversed index.
    pub fn match_to_seq(&self, seq: &[Nucleotide]) -> Vec<PrimerMatch> {
        let mut result = Vec::new();

        // This check prevents spurious small-sequence matches, which may be numerous otherwise.
        if self.sequence.len() < MIN_PRIMER_LEN {
            return result;
        }

        let (matches_fwd, matches_rev) = match_subseq(&self.sequence, seq);

        for range in matches_fwd {
            result.push(PrimerMatch {
                direction: PrimerDirection::Forward,
                range,
            });
        }
        for range in matches_rev {
            result.push(PrimerMatch {
                direction: PrimerDirection::Reverse,
                range,
            });
        }

        result
    }

    /// Formats the indexes, and size of this feature.
    pub fn location_descrip(&self) -> String {
        self.volatile
            .matches
            .iter()
            .map(|match_| format!("{}..{}", match_.range.start, match_.range.end))
            .collect::<Vec<String>>()
            .join("; ")
    }

    /// Automatically select primer length based on quality score.
    pub fn tune(&mut self, ion: &IonConcentrations) {
        match self.volatile.tune_setting {
            TuneSetting::Both(_) => self.tune_both_ends(ion),
            TuneSetting::Disabled => (),
            _ => self.tune_single_end(ion),
        }
    }

    /// Note: In its current form, this assumes only one end is tunable, prior to calling this function.
    fn tune_single_end(&mut self, ion: &IonConcentrations) {
        // todo: Using the seq_input as the only way we store total len feels janky.
        let len_untrimmed = self.volatile.sequence_input.len();

        if len_untrimmed <= MIN_PRIMER_LEN {
            return;
        }

        let mut best_val = 0;
        let mut best_score = 0.;

        let num_vals = len_untrimmed - MIN_PRIMER_LEN;

        // When this function is called, exactly one of these ends must be enabled.
        // let (i) = match &mut self.volatile.tune_setting {
        //     TuneSetting::Only5(v) => v,
        //     TuneSetting::Only3(v) => v,
        //     _ => return,
        // };

        for val in 0..num_vals {
            // We need to have this assignment in the loop to prevent borrow errors.
            let i = match &mut self.volatile.tune_setting {
                TuneSetting::Only5(v) => v,
                TuneSetting::Only3(v) => v,
                _ => return,
            };

            *i = val;
            self.run_calcs(ion);

            if let Some(metrics) = &self.volatile.metrics {
                if metrics.quality_score > best_score {
                    best_val = val;
                    best_score = metrics.quality_score;
                }
            }
        }

        let i = match &mut self.volatile.tune_setting {
            TuneSetting::Only5(v) => v,
            TuneSetting::Only3(v) => v,
            _ => return,
        };
        *i = best_val;
        self.run_calcs(ion);
    }

    fn tune_both_ends(&mut self, ion: &IonConcentrations) {
        // todo: Using the seq_input as the only way we store total len feels janky.
        let len_untrimmed = self.volatile.sequence_input.len();

        // We need the min primer length on both sides of the anchor.
        if len_untrimmed <= MIN_PRIMER_LEN * 2 {
            return;
        }

        let mut best_val = (0, 0);
        let mut best_score = 0.;

        // As for single-ended, we assume this function only runs when both ends are marked tunable.
        let (anchor, _, _) = match self.volatile.tune_setting {
            TuneSetting::Both(v) => v,
            _ => return,
        };

        // We ensure we have the min primer len on either side of the anchor.
        let num_vals_5p = if anchor < MIN_PRIMER_LEN {
            0
        } else {
            anchor - MIN_PRIMER_LEN
        };
        let num_vals_3p = if anchor < MIN_PRIMER_LEN {
            0
        } else {
            if len_untrimmed > anchor {
                let lhs = len_untrimmed - anchor;
                if lhs > MIN_PRIMER_LEN {
                    (len_untrimmed - anchor) - MIN_PRIMER_LEN
                } else {
                    0
                }
            } else {
                0
            }
        };

        // A nested loop: Try all combinations.
        for val5 in 0..num_vals_5p {
            for val3 in 0..num_vals_3p {
                // As for single-ended, we assume this function only runs when both ends are marked tunable.
                let (_, i_5p, i_3p) = match &mut self.volatile.tune_setting {
                    TuneSetting::Both(v) => v,
                    _ => return,
                };

                *i_5p = val5;
                *i_3p = val3;
                self.run_calcs(ion);

                if let Some(metrics) = &self.volatile.metrics {
                    if metrics.quality_score > best_score {
                        best_val = (val5, val3);
                        best_score = metrics.quality_score;
                    }
                }
            }
        }

        let (_, i_5p, i_3p) = match &mut self.volatile.tune_setting {
            TuneSetting::Both(v) => v,
            _ => return,
        };
        *i_5p = best_val.0;
        *i_3p = best_val.1;

        self.run_calcs(ion);
    }

    /// Perform calculations on primer quality and related data. Run this when the sequence changes,
    /// the tuning values change etc.
    ///
    /// This also syncs the active sequence based on the tune settings, and calculates primer weight.
    pub fn run_calcs(&mut self, ion_concentrations: &IonConcentrations) {
        self.volatile.weight = seq_weight(&self.sequence);

        let full_len = self.volatile.sequence_input.len();
        let mut start = 0;
        let mut end = full_len;

        // if let TuneSetting::Enabled(i) = self.volatile.tunable_5p {
        if let Some(i) = self.volatile.tune_setting.val_5p_mut() {
            start = *i;
        }

        // if let TuneSetting::Enabled(i) = self.volatile.tunable_3p {
        if let Some(i) = self.volatile.tune_setting.val_3p_mut() {
            end = if *i > full_len {
                // Prevents an overrun.
                0
            } else {
                full_len - *i
            };
        }

        if start > end || start + 1 > self.volatile.sequence_input.len() {
            start = 0;
            end = full_len
        }

        self.sequence = seq_from_str(&self.volatile.sequence_input[start..end]);
        self.volatile.metrics = self.calc_metrics(ion_concentrations);

        self.volatile.sequence_input[..start].clone_into(&mut self.volatile.seq_removed_5p);
        self.volatile.sequence_input[end..].clone_into(&mut self.volatile.seq_removed_3p);
    }
}

pub fn design_slic_fc_primers(
    seq_vector: &Seq,
    seq_insert: &Seq,
    mut insert_loc: usize,
) -> Option<SlicPrimers> {
    if insert_loc == 0 {
        eprintln!("Error when making SLIC primers: Insert loc is 0");
        return None;
    }
    let seq_len_vector = seq_vector.len();
    let seq_len_insert = seq_insert.len();

    // todo: This likely doesn't handle the off-by-one dynamic.
    insert_loc = insert_loc % seq_len_vector; // Wrap around the origin (circular plasmids)

    let vector_reversed = seq_complement(seq_vector);
    let insert_reversed = seq_complement(seq_insert);

    let insert_loc_reversed = seq_len_vector - insert_loc + 1;

    let (seq_vector_fwd, seq_vector_rev) = {
        let end = (insert_loc + UNTRIMMED_LEN_VECTOR) % seq_len_vector;
        let end_reversed = (insert_loc_reversed + UNTRIMMED_LEN_VECTOR) % seq_len_vector;

        (
            seq_vector[insert_loc - 1..end].to_owned(),
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
            name: "Vector fwd".to_owned(),
            description: Some(
                "SLIC cloning primer, Vector forward. Amplifies the entire vector.".to_owned(),
            ),
            volatile: Default::default(),
        },
        vector_rev: Primer {
            sequence: seq_vector_rev,
            name: "Vector rev".to_owned(),
            description: Some(
                "SLIC cloning primer, Vector reverse. Amplifies the entire vector.".to_owned(),
            ),
            volatile: Default::default(),
        },
        insert_fwd: Primer {
            sequence: seq_insert_fwd,
            name: "Insert fwd".to_owned(),
            description: Some(
                "SLIC cloning primer, Insert forward. Overlaps with the vector.".to_owned(),
            ),
            volatile: Default::default(),
        },
        insert_rev: Primer {
            sequence: seq_insert_rev,
            name: "Insert rev".to_owned(),
            description: Some(
                "SLIC cloning primer, Insert forward. Overlaps with the vector.".to_owned(),
            ),
            volatile: Default::default(),
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
            name: "Amplification fwd".to_owned(),
            description: Some("Amplification primer, forward.".to_owned()),
            volatile: Default::default(),
        },
        rev: Primer {
            sequence: seq_rev,
            name: "Amplification rev".to_owned(),
            description: Some("Amplification primer, reverse.".to_owned()),
            volatile: Default::default(),
        },
    })
}

#[derive(Default, Clone, Encode, Decode)]
pub struct PrimerData {
    /// Editing is handled using this string; we convert the string to our nucleotide sequence as needed.
    /// This includes nts past the tuning offsets.
    pub sequence_input: String,
    pub metrics: Option<PrimerMetrics>,
    // tunable_end: TunableEnd,
    /// These fields control if a given primer end is fixed (Eg marking the start of an insert,
    /// marking the insert point in a vector etc) or if we can tune its length to optimize the primer.
    pub tune_setting: TuneSetting,
    /// These seq_removed fields are redundant with primer and tune settings. We use them to cache
    /// the actual sequences that are removed for display purposes.
    pub seq_removed_5p: String,
    pub seq_removed_3p: String,
    /// todo: Which direction is the range, if the direction is reverse?
    pub matches: Vec<PrimerMatch>,
    /// Primer weight, in Daltons.
    pub weight: f32,
}

impl PrimerData {
    pub fn new(seq: &[Nucleotide]) -> Self {
        let mut result = Self::default();
        result.sequence_input = seq_to_str_lower(seq);
        result
    }
}

// #[derive(Clone, Copy, PartialEq, Debug, Encode, Decode)]
// pub enum TuneSetting {
//     Disabled,
//     /// Inner: Offset index; this marks the distance from the respective ends that the sequence is attentuated to.
//     Enabled(usize),
// }

#[derive(Clone, Copy, PartialEq, Debug, Encode, Decode)]
pub enum TuneSetting {
    /// Inner: Offset index; this marks the distance from the respective ends that the sequence is attentuated to.
    /// The 3' end is the anchor.
    Only5(usize),
    /// The 5' end is the anchor.
    Only3(usize),
    /// There is a central anchor, and both ends are tuanble. We use this for the insert primers of
    /// SLIC and FC.
    Both((usize, usize, usize)), // (Anchor from 5', 5p, 3p)

    Disabled,
    // Enabled(usize),
}

impl Default for TuneSetting {
    fn default() -> Self {
        Self::Disabled
    }
}

impl TuneSetting {
    pub fn toggle_5p(&mut self) {
        *self = match self {
            Self::Only5(_) => Self::Disabled,
            Self::Only3(_) => Self::Both((0, 10, 0)), // todo: 20? Hmm.
            Self::Both(_) => Self::Only3(0),
            Self::Disabled => Self::Only5(0),
        }
    }
    pub fn toggle_3p(&mut self) {
        *self = match self {
            Self::Only3(_) => Self::Disabled,
            Self::Only5(_) => Self::Both((0, 10, 0)), // todo: 20? Hmm.
            Self::Both(_) => Self::Only5(0),
            Self::Disabled => Self::Only3(0),
        }
    }

    /// If the 5p end is tunable, get its valu.
    pub fn val_5p(&self) -> Option<usize> {
        match self {
            Self::Only5(v) => Some(*v),
            Self::Both(v) => Some(v.1),
            _ => None,
        }
    }

    /// If the 5p end is tunable, get its value.
    pub fn val_3p(&self) -> Option<usize> {
        match self {
            Self::Only3(v) => Some(*v),
            Self::Both(v) => Some(v.2),
            _ => None,
        }
    }

    /// If the 5p end is tunable, get its value, mutably.
    pub fn val_5p_mut(&mut self) -> Option<&mut usize> {
        match self {
            Self::Only5(ref mut v) => Some(v),
            Self::Both(ref mut v) => Some(&mut v.1),
            _ => None,
        }
    }

    /// If the 5p end is tunable, get its value, mutably.
    pub fn val_3p_mut(&mut self) -> Option<&mut usize> {
        match self {
            Self::Only3(ref mut v) => Some(v),
            Self::Both(ref mut v) => Some(&mut v.2),
            _ => None,
        }
    }

    pub fn tunable(&self) -> bool {
        !matches!(self, Self::Disabled)
    }
}

/// We run this to generate cloning primers when clicking the button
/// Make sure to do this before inserting the insert into the sequence.
pub fn make_cloning_primers(state: &mut State) {
    let seq_vector = &state.generic[state.active].seq;
    let seq_insert = &state.ui.cloning_insert.seq_insert;

    if let Some(mut primers) =
        design_slic_fc_primers(seq_vector, seq_insert, state.cloning.insert_loc)
    {
        let sequence_input = seq_to_str_lower(&primers.insert_fwd.sequence);

        let insert_fwd_data = PrimerData {
            sequence_input,
            // Both ends are  tunable, since this glues the insert to the vector
            tune_setting: TuneSetting::Both((
                UNTRIMMED_LEN_INSERT,
                DEFAULT_TRIM_AMT,
                DEFAULT_TRIM_AMT,
            )), // todo: TIe the anchor to the const
            ..Default::default()
        };

        let sequence_input = seq_to_str_lower(&primers.insert_rev.sequence);
        let insert_rev_data = PrimerData {
            sequence_input,
            // Both ends are tunable, since this glues the insert to the vector
            tune_setting: TuneSetting::Both((
                UNTRIMMED_LEN_VECTOR,
                DEFAULT_TRIM_AMT,
                DEFAULT_TRIM_AMT,
            )), // todo: QC
            ..Default::default()
        };

        let sequence_input = seq_to_str_lower(&primers.vector_fwd.sequence);
        let vector_fwd_data = PrimerData {
            sequence_input,
            // 5' is non-tunable: This is the insert location.
            tune_setting: TuneSetting::Only3(DEFAULT_TRIM_AMT),
            ..Default::default()
        };

        let sequence_input = seq_to_str_lower(&primers.vector_rev.sequence);
        let vector_rev_data = PrimerData {
            sequence_input,
            // tunable_5p: TuneSetting::Disabled,
            // 3' is non-tunable: This is the insert location.
            tune_setting: TuneSetting::Only3(DEFAULT_TRIM_AMT), // todo: Which one??

            // tunable_3p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
            ..Default::default()
        };

        primers.insert_fwd.volatile = insert_fwd_data;
        primers.insert_rev.volatile = insert_rev_data;
        primers.vector_fwd.volatile = vector_fwd_data;
        primers.vector_rev.volatile = vector_rev_data;

        // Note: If we don't run `run_calcs` before tune here, we get unexpected beavhior culminating
        // in a crash after attempting to tune primers. This is likely related to syncing the tuned-out
        // part of the primers.
        primers
            .insert_fwd
            .run_calcs(&state.ion_concentrations[state.active]);
        primers
            .insert_rev
            .run_calcs(&state.ion_concentrations[state.active]);
        primers
            .vector_fwd
            .run_calcs(&state.ion_concentrations[state.active]);
        primers
            .vector_rev
            .run_calcs(&state.ion_concentrations[state.active]);

        primers
            .insert_fwd
            .tune(&state.ion_concentrations[state.active]);
        primers
            .insert_rev
            .tune(&state.ion_concentrations[state.active]);
        primers
            .vector_fwd
            .tune(&state.ion_concentrations[state.active]);
        primers
            .vector_rev
            .tune(&state.ion_concentrations[state.active]);

        state.generic[state.active].primers.extend([
            primers.insert_fwd,
            primers.insert_rev,
            primers.vector_fwd,
            primers.vector_rev,
        ]);

        state.sync_primer_matches(None);
    }
}

pub fn make_amplification_primers(state: &mut State) {
    if let Some(mut primers) = design_amplification_primers(state.get_seq()) {
        let sequence_input = seq_to_str_lower(&primers.fwd.sequence);

        let primer_fwd_data = PrimerData {
            sequence_input,
            tune_setting: TuneSetting::Only3(DEFAULT_TRIM_AMT),
            // tunable_5p: TuneSetting::Disabled,
            // tunable_3p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
            ..Default::default()
        };

        let sequence_input = seq_to_str_lower(&primers.rev.sequence);
        let primer_rev_data = PrimerData {
            sequence_input,
            tune_setting: TuneSetting::Only3(DEFAULT_TRIM_AMT),
            // tunable_5p: TuneSetting::Disabled,
            // tunable_3p: TuneSetting::Enabled(DEFAULT_TRIM_AMT),
            ..Default::default()
        };

        primers.fwd.volatile = primer_fwd_data;
        primers.rev.volatile = primer_rev_data;

        primers.fwd.tune(&state.ion_concentrations[state.active]);
        primers.rev.tune(&state.ion_concentrations[state.active]);

        state.generic[state.active]
            .primers
            .extend([primers.fwd, primers.rev]);

        state.sync_primer_matches(None); // note: Not requried to run on all primers.
    }
}

#[derive(Clone, Encode, Decode)]
/// Concentrations of common ions in the oglio solution. Affects melting temperature (TM).
/// All values are in milliMolar.
pub struct IonConcentrations {
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
