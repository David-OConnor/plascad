//! Contains code related to identifying toxic proteins.

use na_seq::Seq;

#[derive(Clone, Copy)]
enum Host {
    Ecoli,
    Aav,
}

fn get_toxic_seqs() -> Vec<Seq> {
    vec![]
}
