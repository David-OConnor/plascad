use std::io;

use bincode::{Decode, Encode};

use crate::sequence::Nucleotide::{A, C, G, T};

// Index 0: 5' end.
pub type Seq = Vec<Nucleotide>;

impl Nucleotide {
    pub fn as_str(&self) -> &str {
        match self {
            A => "a",
            T => "t",
            C => "c",
            G => "g",
        }
    }
}

/// A DNA nucleotide.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Encode, Decode)]
#[repr(u8)] // For bio::FASTA compatibility
pub enum Nucleotide {
    A = b'A',
    T = b'T',
    G = b'G',
    C = b'C',
}

impl Nucleotide {
    /// For parsing from bio::FASTA
    pub fn from_u8(val_u8: u8) -> io::Result<Self> {
        match val_u8 {
            b'A' | b'a' => Ok(A),
            b'T' | b't' => Ok(T),
            b'G' | b'g' => Ok(G),
            b'C' | b'c' => Ok(C),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid nucleotide",
            )),
        }
    }
}

#[derive(Clone, Encode, Decode)]
pub struct Feature {
    /// 1-based indexing.
    pub index_range: (usize, usize),
    pub label: String,
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

pub fn seq_to_str(seq: &[Nucleotide]) -> String {
    let mut result = String::new();

    for nt in seq {
        result.push_str(nt.as_str());
    }

    result
}
