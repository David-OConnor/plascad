//! Contains code related to identifying toxic proteins.

use na_seq::Nucleotide;

use crate::protein::proteins_from_seq;

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

fn check_glu_tracts(seq: &[Nucleotide]) -> ToxicStatus {
    // todo: Confirm how you handle reading frames here.
    // let proteins = proteins_from_seq(seq);

    ToxicStatus::Fail
}
