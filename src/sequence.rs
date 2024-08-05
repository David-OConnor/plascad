use std::{collections::HashMap, fmt::Display, io};

use bincode::{Decode, Encode};
use num_enum::TryFromPrimitive;

use crate::{
    gui::navigation::Page,
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

/// A DNA nucleotide. The u8 repr is for use with a compact binary format.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Encode, Decode, TryFromPrimitive)]
#[repr(u8)]
pub enum Nucleotide {
    A = 0,
    T = 1,
    G = 2,
    C = 3,
}

impl Nucleotide {
    /// For parsing from FASTA and SnapGene compatibility.
    pub fn from_u8_letter(val_u8: u8) -> io::Result<Self> {
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

    /// For FASTA and SnapGene compatibility.
    pub fn to_u8_letter(&self) -> u8 {
        match self {
            A => b'A',
            T => b'T',
            G => b'G',
            C => b'C',
        }
    }
}

/// Of the 6 possible reading frames.
#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum ReadingFrame {
    /// Forward, with 0 offset (This pattern applies for all variants)
    Fwd0,
    Fwd1,
    Fwd2,
    Rev0,
    Rev1,
    Rev2,
}

impl ReadingFrame {
    pub fn offset(&self) -> usize {
        match self {
            Self::Fwd0 | Self::Rev0 => 0,
            Self::Fwd1 | Self::Rev1 => 1,
            Self::Fwd2 | Self::Rev2 => 2,
        }
    }
}

impl Default for ReadingFrame {
    fn default() -> Self {
        Self::Fwd0
    }
}

impl Display for ReadingFrame {
    /// For use with selector buttons.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            Self::Fwd0 => "Fwd 0",
            Self::Fwd1 => "Fwd 1",
            Self::Fwd2 => "Fwd 2",
            Self::Rev0 => "Rev 0",
            Self::Rev1 => "Rev 1",
            Self::Rev2 => "Rev 2",
        }
            .to_owned();
        write!(f, "{}", str)
    }
}

#[derive(Debug)]
pub struct ReadingFrameMatch {
    pub range: (usize, usize),
}

#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum FeatureType {
    Generic,
    Gene,
    Ori,
    // RnaPolyBindSite,
    RibosomeBindSite,
    Promoter,
    AntibioticResistance,
    /// Note: This one behaves a bit different from the others; we use it here so we can share the feature
    /// overlay code.
    Primer,
    CodingRegion,
    LongTerminalRepeat,
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
            // Self::RnaPolyBindSite => "RNA poly bind site",
            Self::RibosomeBindSite => "Ribosome bind site",
            Self::Promoter => "Promoter",
            Self::AntibioticResistance => "Antibiotic resistance",
            Self::Primer => "Primer",
            Self::CodingRegion => "Coding region",
            Self::LongTerminalRepeat => "Long term repeat",
        }
            .to_owned()
    }

    pub fn color(&self) -> Color {
        match self {
            Self::Generic => (255, 0, 255),
            Self::Gene => (255, 128, 128),
            Self::Ori => (40, 128, 128),
            // Self::RnaPolyBindSite => (255, 0, 20),
            Self::RibosomeBindSite => (255, 0, 100),
            Self::Promoter => (120, 120, 70),
            Self::AntibioticResistance => (128, 128, 100),
            Self::Primer => (0, 0, 0),             // N/A for now at least.
            Self::CodingRegion => (100, 200, 255), // N/A for now at least.
            Self::LongTerminalRepeat => (150, 200, 255), // N/A for now at least.
        }
    }

    /// Parse from a string; we use this for both SnapGene and GenBank.
    pub fn from_external_str(v: &str) -> Self {
        // todo: Update as required with more
        let v = &v.to_lowercase();

        match v.as_ref() {
            "cds" => Self::CodingRegion,
            "rbs" => Self::RibosomeBindSite,
            "rep_origin" => Self::Ori,
            "promoter" => Self::Promoter,
            "primer_bind" => Self::Primer, // todo: This is a bit awk; genbank.
            "ltr" => Self::LongTerminalRepeat,
            "misc_feature" => Self::Generic,
            _ => Self::Generic,
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
    pub notes: HashMap<String, String>,
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
    let stop_codons = vec![[T, A, A], [T, A, G], [T, G, A]];

    const START_CODON: [Nucleotide; 3] = [A, T, G];

    let mut result = Vec::new();

    result
}
