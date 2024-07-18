//! This module contains code related to primer (oglionucleotide) design and QC.

use bincode::{Decode, Encode};

use crate::{
    Nucleotide::{C, G},
    Seq,
    util::{map_linear, seq_complement},
};
use crate::util::seq_from_str;

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

#[derive(Clone, Copy, Encode, Decode)]
pub enum PrimerDirection {
    Forward,
    Reverse,
}

/// Metrics related to primer quality.
#[derive(Clone, Debug, Default, Encode, Decode)]
pub struct PrimerMetrics {
    /// C
    pub melting_temp: f32,
    /// 0. to 1.
    pub gc_portion: f32,
    /// How many G and C nts are in the last 5 (3' end) nts of the sequence.
    pub gc_3p_count: u8,
    pub complexity: f32,
    pub self_end_dimer: u8,
    pub repeats: u8,
    pub tm_score: f32,
    pub gc_score: f32,
    pub gc_3p_score: f32,
    pub complexity_score: f32,
    pub dimer_score: f32,
    /// https://www.benchling.com/primer-design-for-pcr
    /// "Avoid runs of four or more of a single base (e.g., ACCCCC), or four or more dinucleotide
    /// repeats (e.g., ATATATATAT) as they will cause mispriming."
    pub repeats_score: f32,
    /// This is a weighted overall score.
    pub quality_score: f32,
}

impl PrimerMetrics {
    /// Return a quality score, on a scale from 0 to 1.
    pub fn update_scores(&mut self) {
        const GC_TARGET: f32 = 0.5;

        const WEIGHT_TM: f32 = 1.;
        const WEIGHT_GC: f32 = 1.;
        const WEIGHT_STAB: f32 = 1.;
        const WEIGHT_COMPLEXITY: f32 = 1.;
        const WEIGHT_DIMER: f32 = 1.;
        const WEIGHT_REPEATS: f32 = 1.;

        // todo: Instead of closeness to 59, should it be >54??
        // Also: 50-60C. And within 5C of the complement primer.
        self.tm_score = map_linear((self.melting_temp - TM_TARGET).abs(), (0., 12.), (1., 0.));
        self.tm_score = self.tm_score.clamp(0., 1.);

        // This is currently a linear map, between 0 and 1.
        self.gc_score = 1. - (self.gc_portion - GC_TARGET).abs() * 2.;

        self.gc_3p_score = match self.gc_3p_count {
            0 | 1 => 1.,
            2 => 0.9,
            3 => 0.4,
            4 => 0.1,
            5 => 0.,
            _ => unreachable!(),
        };

        self.repeats_score = match self.repeats {
            0 => 1.,
            1 => 0.8,
            2 => 0.4,
            3 => 0.1,
            _ => 0.,
        };

        self.quality_score = ((WEIGHT_TM * self.tm_score
            + WEIGHT_GC * self.gc_score
            + WEIGHT_STAB * self.gc_3p_score
            + WEIGHT_COMPLEXITY * self.complexity_score
            + WEIGHT_DIMER * self.dimer_score)
            + WEIGHT_REPEATS * self.repeats_score) / 6.
    }
}

#[derive(Default, Encode, Decode)]
pub struct Primer {
    pub sequence: Seq,
}

impl Primer {
    /// Calculate melting temperature (TM), in C.
    ///
    /// [Calc data from AmplifX](https://inp.univ-amu.fr/en/amplifx-manage-test-and-design-your-primers-for-pcr):
    /// "TM = 81.5 +16.6 X log10([Na+]+[K+])+0.41 x(%GC) - 675/N (with default values: [ Na+]+[K+]=0.05 (50mM))
    ///
    /// "the most precise... "bases stacking method” : TM=dH/(dS+0.368xNxln[Na+]+ R
    /// x ln[Primer]/4) with R=1.987 and the different dH and dS taken in [Santalucia J PNAS 95 pp1460-1465 (1998)]"
    ///
    /// https://primerexplorer.jp/e/v4_manual/03.html (Different formula)
    ///
    /// https://www.rosalind.bio/en/knowledge/what-formula-is-used-to-calculate-tm
    /// Come up with something better A/R
    /// "ASSUMPTIONS:
    /// Equations above assume that the annealing occurs under the standard conditions of 50 nM
    /// primer, 50 mM Na+, and pH 7.0."
    pub fn calc_tm(&self) -> f32 {
        // https://www.rosalind.bio/en/knowledge/what-formula-is-used-to-calculate-tm
        // Come up with something better A/R
        // "ASSUMPTIONS:
        // Equations above assume that the annealing occurs under the standard conditions of 50 nM
        // primer, 50 mM Na+, and pH 7.0."

        // let mut num_a = 0;
        // let mut num_c = 0;
        // let mut num_t = 0;
        // let mut num_g = 0;
        //
        // for nt in &self.sequence {
        //     match nt {
        //         A => num_a += 1,
        //         T => num_t += 1,
        //         C => num_c += 1,
        //         G => num_g += 1,
        //     }
        // }
        //
        // if self.sequence.len() < 14 {
        //     ((num_a + num_t) * 2 + (num_g + num_c) * 4) as f32
        // } else {
        //     64.9 + 41. * (num_g as f32 + num_c as f32 - 16.4) / (num_a + num_t + num_g + num_c) as f32
        // }

        const NAP_K_P: f32 = 0.05;

        // TM = 81.5 +16.6 X log10([Na+]+[K+])+0.41 x(%GC) - 675/N
        81.5 + 16.6 * NAP_K_P.log10() + 0.41 * self.calc_gc() * 100.
            - 675. / (self.sequence.len() as f32)
    }

    /// Calculate GC portion, on a scale of 0 to 1.
    pub fn calc_gc(&self) -> f32 {
        let mut num_gc = 0;
        for nt in &self.sequence {
            if *nt == C || *nt == G {
                num_gc += 1;
            }
        }

        num_gc as f32 / self.sequence.len() as f32
    }

    /// This is a metric known as 3' end stability. Return the number of Gs and Cs in the last 5 bases.
    /// A target is no more than 2.
    pub fn count_3p_g_c(&self) -> u8 {
        let mut result = 0;

        let len = self.sequence.len();
        let last_5 = if len > 5 {
            self.sequence.split_at(len - 5).1
        } else {
            &self.sequence[..]
        };

        for nt in last_5 {
            if *nt == C || *nt == G {
                result += 1
            }
        }

        result
    }

    /// Calculate complexity. (Details on what this means here)
    pub fn calc_complexity(&self) -> f32 {
        0.
    }

    /// Calculate self end dimer score. (Details on what this means here)
    ///
    /// http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
    /// "A primer self-dimer is formed by intermolecular interactions between the two (same sense)
    /// primers, where the primer is homologous to itself. Generally a large amount of primers are
    /// used in PCR compared to the amount of target gene. When primers form intermolecular dimers
    /// much more readily than hybridizing to target DNA, they reduce the product yield. Optimally a 3'
    /// end self dimer with a ΔG of -5 kcal/mol and an internal self dimer with a ΔG of -6 kcal/mol is
    /// tolerated generally."
    pub fn calc_self_end_dimer(&self) -> u8 {
        0
    }

    /// Calculate how many single or double nucleotide sequences exist that are of len 4 or more of the
    /// same nt or nt pair respectively. It currently doesn't differentiate between 4, and more.
    pub fn calc_repeats(&self) -> u8 {
        if self.sequence.len() < 4 {
            return 0;
        }
        
        
        let mut result = 0;

        // Check for single nt repeats.
        let mut prev_nt = self.sequence[0];
        let mut single_repeat_len = 0;
        
        for nt in &self.sequence {
            println!("NT: {:?}, prev_nt: {:?}", nt, prev_nt);
            if *nt == prev_nt {
                single_repeat_len += 1;
                
                if single_repeat_len > 2 { // todo: Why off by one?
                    result += 1;
                    single_repeat_len = 0;
                }
            } else {
                single_repeat_len = 0;
            }
            prev_nt = *nt;
        }

        let mut prev_nts = (self.sequence[0], self.sequence[1]);
        let mut double_repeat_len = 0;

        // todo: Come back to adn implement this.
        // for nt in &self.sequence {
        //     if nt == prev_nt {
        //         prev_nt = nt;
        //         double_repeat_len += 1;
        //
        //         if double_repeat_len >= 4 {
        //             result += 1;
        //             double_repeat_len = 0;
        //         }
        //     }
        // }
        
        result
    }

    /// Calculate all primer metrics.
    /// todo: methods on Metrics instead?
    pub fn calc_metrics(&self) -> Option<PrimerMetrics> {
        if self.sequence.len() < MIN_PRIMER_LEN {
            return None;
        }

        let mut result = PrimerMetrics {
            melting_temp: self.calc_tm(),
            gc_portion: self.calc_gc(),
            gc_3p_count: self.count_3p_g_c(),
            complexity: self.calc_complexity(),
            self_end_dimer: self.calc_self_end_dimer(),
            repeats: self.calc_repeats(),
            ..Default::default()
        };
        result.update_scores();

        Some(result)
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
        // let insert_only_seq_rev: Vec<Nucleotide> = Vec::new();

        let mut fwd = seq_vector_fwd.clone();
        fwd.extend(insert_only_seq_fwd);

        let mut rev = seq_vector_fwd.clone();
        rev.extend(insert_only_seq_rev);

        (fwd, rev)
    };

    let mut result = SlicPrimers {
        vector_fwd: Primer {
            sequence: seq_vector_fwd,
        },
        vector_rev: Primer {
            sequence: seq_vector_rev,
        },
        insert_fwd: Primer {
            sequence: seq_insert_fwd,
        },
        insert_rev: Primer {
            sequence: seq_insert_rev,
        },
    };

    // todo: Optimize.

    Some(result)
}

// todo: Use this A/R, called from the UI page.
pub fn design_amplification_primers(seq: &Seq) -> Option<AmplificationPrimers> {
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

    let mut result = AmplificationPrimers {
        fwd: Primer { sequence: seq_fwd },
        rev: Primer { sequence: seq_rev },
    };

    // todo: Optimize.

    Some(result)
}

#[derive(Default, Encode, Decode)]
pub struct PrimerData {
    /// This primer excludes nts past the tuning offsets.
    pub primer: Primer,
    /// Editing is handled using this string; we convert the string to our nucleotide sequence as needed.
    /// This includes nts past the tuning offsets.
    pub sequence_input: String,
    pub description: String,
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
}

impl PrimerData {
    /// Perform calculations on primer quality and related data. Run this when the sequence changes,
    /// the tuning values change etc.
    pub fn run_calcs(&mut self) {
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
        self.metrics = self.primer.calc_metrics();

        self.sequence_input[..start].clone_into(&mut self.seq_removed_5p);
        self.sequence_input[end..].clone_into(&mut self.seq_removed_3p);
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