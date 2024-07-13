//! This module contains code related to primer design and QC

use bio::io::fasta::Sequence;
use crate::{Nucleotide, Seq};

struct SlicPrimers {
    pub vector_fwd: Seq,
    pub vector_rev: Seq,
    pub insert_fwd: Seq,
    pub insert_rev: Seq,
}

/// Metrics related to primer quality.
struct PrimerMetrics {
    /// C
    tm: f32,
    /// 0. to 1.
    gc_portion: f32,
    end_stability: u8,
    complexity: f32,
    self_end_dimer: u8,
}

impl PrimerMetrics {
    pub fn score(&self) {
        const GC_TARGET: f32 = 0.5;

        // This is currently a linear map, between 0 and 1.
        let score_gc = 1. - (self.gc_portion - GC_TARGET).abs() * 2.;
    }
}

struct Primer {
    nucleotides: Seq
}

impl Primer {
    /// Calculate melting temperature (TM), in C.
    pub fn calc_tm(&self) -> f32 {

        // https://www.rosalind.bio/en/knowledge/what-formula-is-used-to-calculate-tm
        // Come up with something better A/R
        // "ASSUMPTIONS:
        // Equations above assume that the annealing occurs under the standard conditions of 50 nM
        // primer, 50 mM Na+, and pH 7.0."

        let mut num_a = 0;
        let mut num_c = 0;
        let mut num_t = 0;
        let mut num_g = 0;

        for nt in &self.nucleotides {
            match nt {
                Nucleotide::A => num_a += 1,
                Nucleotide::T => num_t += 1,
                Nucleotide::C => num_c += 1,
                Nucleotide::G => num_g += 1,
            }
        }

        if self.nucleotides.len() < 14 {
            ((num_a + num_t) * 2 + (num_g + num_c) * 4) as f32
        } else {
            64.9 + 41. * (num_g as f32 + num_c as f32 - 16.4) / (num_a + num_t + num_g + num_c) as f32
        }
    }

    /// Calculate GC portion, on a scale of 0 to 1.
    pub fn calc_gc(&self) -> f32 {
        let mut num_gc = 0;
        for nt in &self.nucleotides {
            if nt == *Nucleotide::C || nt == *Nucleotide::G {
                num_gc += 1;
            }
        }

        num_gc as f32 / self.nucleotides.len() as f32
    }

    /// Calculate 3' end stability. (Details on what this means here)
    pub fn calc_end_stability(&self) -> f32 {

    }

    /// Calculate complexity. (Details on what this means here)
    pub fn calc_complexity(&self) -> f32 {
        let mut num_gc = 0;
        for nt in &self.nucleotides {
            if nt == *Nucleotide::C || nt == *Nucleotide::G {
                num_gc += 1;
            }
        }

        num_gc as f32 / self.nucleotides.len() as f32
    }

    /// Calculate self end dimer score. (Details on what this means here)
    pub fn calc_self_send_dimer(&self) -> f32 {
    }

    /// Calculate all primer metrics.
    pub fn calc_metrics(&self) -> PrimerMetrics {
        PrimerMetrics {
            tm: self.calc_tm(),
            gc_portion: self.calc_gc(),
            end_stability: self.calc_end_stability(),
            complexity: self.calc_complexity(),
            self_end_dimer: self.calc_self_end_dimer(),
        }
    }
}

pub fn design_slic_fc_primers(seq_vector: &Seq, insert_i: usize, seq_insert: &Seq) -> SlicPrimers {

}