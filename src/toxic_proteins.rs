//! Contains code related to identifying toxic proteins.

use na_seq::{
    amino_acids::{AminoAcid, CodingResult},
    Nucleotide,
};

use crate::protein::proteins_from_seq;

const GLN_TRACT_SCAN_FRAME: usize = 35;
const GLN_MAX_PORTION: f32 = 0.75;

#[derive(Clone, Copy, PartialEq)]
enum ToxicStatus {
    Pass,
    Fail,
}

// #[derive(Clone, Copy)]
// enum Host {
//     Ecoli,
//     Aav,
// }

/// Identify tracts, in all reading frames, of long sequences of Glutamine; this may lead to aggregation,
/// and therefor toxicity. Sequences of 35 Gln in a row are a good indicator of this.
fn check_glu_tracts(seq: &[Nucleotide]) -> ToxicStatus {
    let len = seq.len();

    let mut seq_rev = seq.to_vec();
    seq_rev.reverse();

    // For each reading frame, create a Vec of amino acids.
    for seq_ in &[seq, &seq_rev] {
        for rf_offset in 0..2 {
            let mut aa_seq = Vec::new();

            for i_ in 0..len / 3 {
                let i = i_ * 3 + rf_offset; // The ORF-modified sequence index.

                if i + 3 >= len {
                    break;
                }

                let nts = &seq_[i..i + 3];

                // let mut matched = false;
                if let CodingResult::AminoAcid(aa) = AminoAcid::from_codons(nts.try_into().unwrap())
                {
                    // Note: We are ignoring stop codons here.
                    aa_seq.push(aa)
                }
            }

            // Now  that we have a set of AA for this reading frame, scan for Glu tracts.
            // If we find any, return the failure state.

            for start_i in 0..aa_seq.len() - GLN_TRACT_SCAN_FRAME {
                let mut num_gln = 0;
                for i in start_i..start_i + aa_seq.len() {
                    if aa_seq[i] == AminoAcid::Gln {
                        num_gln += 1;
                    }
                }

                if num_gln as f32 / GLN_TRACT_SCAN_FRAME as f32 > GLN_MAX_PORTION {
                    return ToxicStatus::Fail;
                }
            }
        }
    }

    ToxicStatus::Pass
}
