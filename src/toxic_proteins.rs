//! Contains code related to identifying toxic proteins.

use seq::Seq;

#[derive(Clone, Copy)]
enum Host {
    Ecoli,
    Aav,
}

fn get_toxic_seqs() -> Vec<Seq> {
    vec![]
}
