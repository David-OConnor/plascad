//! This module contains code related to primer design and QC

// use bio::io::fasta::Sequence;

use crate::{
    util::map_linear,
    Nucleotide::{self, A, C, G, T},
    Seq,
};

pub const MIN_PRIMER_LEN: usize = 5;

struct SlicPrimers {
    pub vector_fwd: Seq,
    pub vector_rev: Seq,
    pub insert_fwd: Seq,
    pub insert_rev: Seq,
}

/// Metrics related to primer quality.
#[derive(Clone, Debug)]
pub struct PrimerMetrics {
    /// C
    pub melting_temp: f32,
    /// 0. to 1.
    pub gc_portion: f32,
    pub gc_3p_count: u8,
    pub complexity: f32,
    pub self_end_dimer: u8,
}

impl PrimerMetrics {
    /// Return a quality score, on a scale from 0 to 1.
    pub fn quality(&self) -> f32 {
        const TM_TARGET: f32 = 59.;
        const GC_TARGET: f32 = 0.5;

        const WEIGHT_TM: f32 = 1.;
        const WEIGHT_GC: f32 = 1.;
        const WEIGHT_STAB: f32 = 1.;
        const WEIGHT_COMPLEXITY: f32 = 1.;
        const WEIGHT_DIMER: f32 = 1.;

        let score_tm = map_linear((self.melting_temp - TM_TARGET).abs(), (0., 12.), (1., 0.));
        let score_tm = score_tm.clamp(0., 1.);

        // This is currently a linear map, between 0 and 1.
        let score_gc = 1. - (self.gc_portion - GC_TARGET).abs() * 2.;

        let score_3p_end_stab = match self.gc_3p_count {
            0 | 1 => 1.,
            2 => 0.9,
            3 => 0.4,
            4 => 0.1,
            5 => 0.,
            _ => unreachable!(),
        };

        let score_complexity = 1.;

        let score_dimer = 1.;

        WEIGHT_TM * score_tm
            + WEIGHT_GC
            + score_gc
            + WEIGHT_STAB * score_3p_end_stab
            + WEIGHT_COMPLEXITY * score_complexity
            + WEIGHT_DIMER * score_dimer
    }
}

#[derive(Default)]
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

    /// Calculate all primer metrics.
    pub fn calc_metrics(&self) -> Option<PrimerMetrics> {
        if self.sequence.len() < MIN_PRIMER_LEN {
            return None;
        }

        Some(PrimerMetrics {
            melting_temp: self.calc_tm(),
            gc_portion: self.calc_gc(),
            gc_3p_count: self.count_3p_g_c(),
            complexity: self.calc_complexity(),
            self_end_dimer: self.calc_self_end_dimer(),
        })
    }
}

pub fn design_slic_fc_primers(seq_vector: &Seq, insert_i: usize, seq_insert: &Seq) -> SlicPrimers {
    SlicPrimers {
        vector_fwd: Vec::new(),
        vector_rev: Vec::new(),
        insert_fwd: Vec::new(),
        insert_rev: Vec::new(),
    }
}
