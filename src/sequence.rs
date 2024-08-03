use std::io;

use bincode::{Decode, Encode};

use crate::{
    primer::PrimerDirection,
    sequence::Nucleotide::{A, C, G, T},
    Color, StateUi,
};

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

#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum FeatureType {
    Generic,
    Gene,
    Ori,
    RnaPolyBindSite,
    AntibioticResistance,
    /// Note: This one behaves a bit different from the others; we use it here so we can share the feature
    /// overlay code.
    Primer,
}

impl Default for FeatureType {
    fn default() -> Self {
        Self::Generic
    }
}

impl FeatureType {
    pub fn to_string(&self) -> String {
        match self {
            Self::Generic => "Generic",
            Self::Gene => "Gene",
            Self::Ori => "Origin of replication",
            Self::RnaPolyBindSite => "RNA polymerase bind site",
            Self::AntibioticResistance => "Antibiotic resistance",
            Self::Primer => "Primer",
        }
        .to_owned()
    }

    pub fn color(&self) -> Color {
        match self {
            Self::Generic => (255, 0, 255),
            Self::Gene => (255, 128, 128),
            Self::Ori => (40, 128, 128),
            Self::RnaPolyBindSite => (255, 0, 20),
            Self::AntibioticResistance => (128, 128, 100),
            Self::Primer => (0, 0, 0), // N/A for now at least.
        }
    }
}

#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum FeatureDirection {
    None,
    Forward,
    Reverse,
}

impl From<PrimerDirection> for FeatureDirection {
    fn from(value: PrimerDirection) -> Self {
        match value {
            PrimerDirection::Forward => Self::Forward,
            PrimerDirection::Reverse => Self::Reverse,
        }
    }
}

impl Default for FeatureDirection {
    fn default() -> Self {
        Self::None
    }
}

impl FeatureDirection {
    pub fn to_string(&self) -> String {
        match self {
            Self::None => "None",
            Self::Forward => "Forward",
            Self::Reverse => "Reverse",
        }
        .to_owned()
    }
}

#[derive(Clone, Encode, Decode)]
pub struct Feature {
    /// 1-based indexing.
    pub index_range: (usize, usize),
    pub feature_type: FeatureType,
    pub direction: FeatureDirection,
    pub label: String,
    /// By default, we display features using featuretype-specific color. Allow the user
    /// to override this.
    pub color_override: Option<Color>,
}

#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum SeqTopology {
    Linear,
    Circular,
}

impl Default for SeqTopology {
    fn default() -> Self {
        Self::Circular
    }
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

// todo: This may not be feasible without aligning reading frames. Come back to this.
/// Automatically generate transient features for start and stop codons.
pub fn _start_stop_codons(seq: &[Nucleotide]) -> Vec<Feature> {
    let start_codons = vec![[T, A, A], [T, A, G], [T, G, A]];

    const START_CODON: [Nucleotide; 3] = [A, T, G];

    let mut result = Vec::new();

    result
}
