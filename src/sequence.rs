use bincode::{Decode, Encode};

use crate::sequence::Nucleotide::{A, C, G, T};

impl Nucleotide {
    pub fn as_str(&self) -> &str {
        match self {
            Self::A => "a",
            Self::T => "t",
            Self::C => "c",
            Self::G => "g",
        }
    }
}

/// A DNA nucleotide.
/// todo: RNA A/R
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Encode, Decode)]
pub enum Nucleotide {
    A,
    T,
    G,
    C,
}

// Index 0: 5' end.
pub type Seq = Vec<Nucleotide>;

#[derive(Encode, Decode)]
pub struct Feature {
    pub index_range: (usize, usize),
    pub name: String,
    pub color: (u8, u8, u8),
}

/// Reverse direction, and swap C for G, A for T.
pub fn seq_complement(seq: &[Nucleotide]) -> Seq {
    let mut result = seq.to_vec();
    result.reverse();

    for nt in &mut result {
        *nt = match *nt {
            A => T,
            T => A,
            C => G,
            G => C,
        };
    }

    result
}

pub fn seq_from_str(str: &str) -> Seq {
    let mut result = Vec::new();

    for char in str.to_lowercase().chars() {
        match char {
            'a' => result.push(A),
            't' => result.push(T),
            'c' => result.push(C),
            'g' => result.push(G),
            _ => (),
        };
    }

    result
}

pub fn make_seq_str(seq: &[Nucleotide]) -> String {
    let mut result = String::new();

    for nt in seq {
        result.push_str(nt.as_str());
    }

    result
}
