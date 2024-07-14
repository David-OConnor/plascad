//! Contains code related to identifying toxic proteins.

// use bio::bio_types::sequence::Sequence;

use crate::Seq;

enum AminoAcid {
    A,
}

#[derive(Clone, Copy)]
enum Host {
    Ecoli,
    Aav,
}

fn get_toxic_seqs() -> Vec<Seq> {
    vec![]
}
