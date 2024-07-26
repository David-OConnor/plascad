//! This module handles assessing various primer metrics, such as GC concentration, and repeats.

use bincode::{Decode, Encode};

use crate::{
    melting_temp_calcs,
    primer::{calc_gc, Primer, MIN_PRIMER_LEN},
    util::{map_linear, remove_duplicates},
    IonConcentrations,
    Nucleotide::{C, G},
    TM_TARGET,
};

/// Metrics related to primer quality.
#[derive(Clone, Debug, Default, Encode, Decode)]
pub struct PrimerMetrics {
    /// C
    pub melting_temp: f32,
    /// 0. to 1.
    pub gc_portion: f32,
    /// How many G and C nts are in the last 5 (3' end) nts of the sequence.
    pub gc_3p_count: u8,
    pub self_end_dimer: u8,
    pub repeats: u8,
    pub tm_score: f32,
    pub gc_score: f32,
    pub gc_3p_score: f32,
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
        // const WEIGHT_COMPLEXITY: f32 = 1.;
        const WEIGHT_DIMER: f32 = 1.;
        const WEIGHT_REPEATS: f32 = 1.;

        // todo: Instead of closeness to 59, should it be >54??
        // Also: 50-60C. And within 5C of the complement primer.
        self.tm_score = map_linear((self.melting_temp - TM_TARGET).abs(), (0., 12.), (1., 0.));
        self.tm_score = self.tm_score.clamp(0., 1.);

        // This is currently a linear map, between 0 and 1.
        self.gc_score = 1. - (self.gc_portion - GC_TARGET).abs() * 2.;

        // Sources differ on if 4 is an ok value. AmplifX calls it "good"; [the DNA universe](https://the-dna-universe.com/2022/09/05/primer-design-guide-the-top-5-factors-to-consider-for-optimum-performance/)
        // considers it to be bad.
        self.gc_3p_score = match self.gc_3p_count {
            0 | 1 => 0.,
            2 => 1.,
            3 => 1.,
            4 => 0.5,
            5 => 0.,
            _ => unreachable!(),
        };

        self.repeats_score = match self.repeats {
            0 => 1.,
            1 => 0.8,
            2 => 0.6,
            3 => 0.3,
            4 => 0.2,
            _ => 0.,
        };

        self.quality_score = ((WEIGHT_TM * self.tm_score
            + WEIGHT_GC * self.gc_score
            + WEIGHT_STAB * self.gc_3p_score
            // + WEIGHT_COMPLEXITY * self.complexity_score
            + WEIGHT_DIMER * self.dimer_score)
            + WEIGHT_REPEATS * self.repeats_score)
            / 5.
    }
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
    ///
    /// We use the "bases stacking method" as defined hte Amplifx guide above.
    ///
    /// See [The BioPython MeltingTemp module](https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/MeltingTemp.py)
    pub fn calc_tm(&self, ion_concentrations: &IonConcentrations) -> f32 {
        // const NAP_K_P: f32 = 0.05;
        //
        // // TM = 81.5 +16.6 X log10([Na+]+[K+])+0.41 x(%GC) - 675/N
        // // 81.5 + 16.6 * NAP_K_P.log10() + 0.41 * self.calc_gc() * 100.
        // //     - 675. / (self.sequence.len() as f32)
        //
        // const dH: f32 = 1. // todo
        // const dS: f32 = 1. // todo
        // const R: f32 =  1.987; // Universal gas constant (Cal/C * Mol)
        // let N = self.sequence.len();
        //
        // dH / (dS + 0.368 * N * NAP.ln() + R * primer.ln() / 4.)
        melting_temp_calcs::calc_tm(&self.sequence, ion_concentrations).unwrap_or(0.)
    }

    /// This is a metric known as 3' end stability. Return the number of Gs and Cs in the last 5 bases.
    /// 2-4 is ideal. (Other sources, 2-3, with 4 being too many). 1 or 5 is considered suboptimal.
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
    /// Also includes repeats of a set of 3 nucleotides anywhere in the seq.
    pub fn calc_repeats(&self) -> u8 {
        if self.sequence.len() < 4 {
            return 0;
        }

        let mut result = 0;

        // Check for single nt repeats.
        let mut prev_nt = self.sequence[0];
        let mut single_repeat_len = 0;

        for nt in &self.sequence {
            if *nt == prev_nt {
                single_repeat_len += 1;

                if single_repeat_len > 2 {
                    // todo: Why off by one?
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

        let mut triplets = Vec::new();
        for (i, nt) in self.sequence.iter().enumerate() {
            if i == self.sequence.len() - 2 {
                break;
            }
            triplets.push((i, nt, self.sequence[i + 1], self.sequence[i + 2]));
        }

        let mut triplet_repeat_seqs = Vec::new();

        for (i, nt) in self.sequence.iter().enumerate() {
            if i == self.sequence.len() - 2 {
                break;
            }

            let triplet_this = (nt, self.sequence[i + 1], self.sequence[i + 2]);

            for triplet_other in &triplets {
                if triplet_this == (triplet_other.1, triplet_other.2, triplet_other.3)
                    && i != triplet_other.0
                {
                    // Dount count each additional nt match beyond 3 as a new repeat
                    if i >= 3
                        && triplet_other.0 >= 3
                        && self.sequence[i - 1] == self.sequence[triplet_other.0 - 1]
                    {
                        continue;
                    }

                    triplet_repeat_seqs.push(triplet_this);

                    // println!("Triplet rep: i this: {}, i_other: {:?} seq: {:?}", i, triplet_other.0, (triplet_other.1, triplet_other.2, triplet_other.3));
                    // triple_nt_repeats += 1;
                }
            }
        }

        triplet_repeat_seqs = remove_duplicates(triplet_repeat_seqs);
        // println!("TRS: {:?}, ", triplet_repeat_seqs);
        result += triplet_repeat_seqs.len() as u8;

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
    pub fn calc_metrics(&self, ion_concentrations: &IonConcentrations) -> Option<PrimerMetrics> {
        if self.sequence.len() < MIN_PRIMER_LEN {
            return None;
        }

        let mut result = PrimerMetrics {
            melting_temp: self.calc_tm(ion_concentrations),
            gc_portion: calc_gc(&self.sequence),
            gc_3p_count: self.count_3p_g_c(),
            // complexity: self.calc_complexity(),
            self_end_dimer: self.calc_self_end_dimer(),
            repeats: self.calc_repeats(),
            ..Default::default()
        };
        result.update_scores();

        Some(result)
    }
}
