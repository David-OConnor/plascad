use std::{collections::HashMap, fmt::Display, io};

use bincode::{Decode, Encode};
use num_enum::TryFromPrimitive;

use crate::{
    primer::PrimerDirection,
    sequence::Nucleotide::{A, C, G, T},
    Color,
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
#[derive(Clone, Copy, PartialEq, Debug, Encode, Decode)]
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
    pub frame: ReadingFrame,
    /// 1-based indexing.
    /// Indices are respective to the non-complementary seq, for both forward and reverse reading frames.
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
    /// Ie, a gene.
    CodingRegion,
    LongTerminalRepeat,
    /// We don't draw these on the map or sequence views; found in GenBank formats (at least), these
    /// are the range of the entire sequence.
    Source,
    Exon,
    Transcript,
}

impl Default for FeatureType {
    fn default() -> Self {
        Self::Generic
    }
}

impl FeatureType {
    /// For displaying in the UI
    pub fn to_string(&self) -> String {
        match self {
            Self::Generic => "Generic",
            Self::Gene => "Gene",
            Self::Ori => "Origin of replication",
            Self::RibosomeBindSite => "Ribosome bind site",
            Self::Promoter => "Promoter",
            Self::AntibioticResistance => "Antibiotic resistance",
            Self::Primer => "Primer",
            Self::CodingRegion => "Coding region",
            Self::LongTerminalRepeat => "Long term repeat",
            Self::Source => "Source",
            Self::Exon => "Exon",
            Self::Transcript => "Transcript",
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
            Self::Source => (120, 70, 120),
            Self::Exon => (255, 255, 180),
            Self::Transcript => (180, 255, 180),
        }
    }

    /// Parse from a string; we use this for both SnapGene and GenBank.
    pub fn from_external_str(v: &str) -> Self {
        // todo: Update as required with more
        let v = &v.to_lowercase();

        match v.as_ref() {
            "cds" => Self::CodingRegion,
            "gene" => Self::Gene,
            "rbs" => Self::RibosomeBindSite,
            "rep_origin" => Self::Ori,
            "promoter" => Self::Promoter,
            "primer_bind" => Self::Primer, // todo: This is a bit awk; genbank.
            "ltr" => Self::LongTerminalRepeat,
            "misc_feature" => Self::Generic,
            "source" => Self::Source,
            "exon" => Self::Exon,
            "transcript" => Self::Transcript,
            _ => Self::Generic,
        }
    }

    /// Create a string for use with SnapGene and GenBank formats.
    pub fn to_external_str(&self) -> String {
        // todo: Update as required with more
        match self {
            Self::Generic => "misc_feature",
            Self::Gene => "gene",
            Self::Ori => "rep_origin",
            Self::RibosomeBindSite => "rbs",
            Self::Promoter => "promoter",
            Self::AntibioticResistance => "antibiotic resistance", // todo
            Self::Primer => "primer_bind",
            Self::CodingRegion => "cds",
            Self::LongTerminalRepeat => "ltr",
            Self::Source => "source",
            Self::Exon => "exon",
            Self::Transcript => "transcript",
        }
        .to_string()
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
    /// 1-based indexing, inclusive. (Note: Could also use the builtin RangeInclusive.)
    pub index_range: (usize, usize),
    pub feature_type: FeatureType,
    pub direction: FeatureDirection,
    pub label: String,
    /// By default, we display features using featuretype-specific color. Allow the user
    /// to override this.
    pub color_override: Option<Color>,
    pub notes: HashMap<String, String>,
}

impl Feature {
    /// Get the color to draw; type color, unless overridden.
    pub fn color(&self) -> Color {
        match self.color_override {
            Some(c) => c,
            None => self.feature_type.color()
        }
    }
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

/// Convert a nucleotide sequence to string.
pub fn seq_to_str(seq: &[Nucleotide]) -> String {
    let mut result = String::new();

    for nt in seq {
        result.push_str(nt.as_str());
    }

    result
}

/// Find coding regions in a sequence, given a reading frame.
pub fn find_orf_matches(seq: &[Nucleotide], orf: ReadingFrame) -> Vec<ReadingFrameMatch> {
    const START_CODON: [Nucleotide; 3] = [A, T, G];
    let stop_codons = [[T, A, A], [T, A, G], [T, G, A]];

    let mut result = Vec::new();

    let mut offset = orf.offset();

    let seq_len_full = seq.len();

    let seq = &match orf {
        ReadingFrame::Fwd0 | ReadingFrame::Fwd1 | ReadingFrame::Fwd2 => seq.to_vec(),
        _ => seq_complement(seq),
    }[offset..];

    let len = seq.len();

    let mut frame_open = None; // Inner: Start index.

    for i_ in 0..len / 3 {
        let i = i_ * 3; // The actual sequence index.

        let nts = &seq[i..i + 3];

        if frame_open.is_none() && nts == START_CODON {
            frame_open = Some(i);

            // println!("Frame open: {:?}", i);
        } else if frame_open.is_some() && stop_codons.contains(nts.try_into().unwrap()) {
            // println!("Frame cl: {:?}", i);

            // + 1 for our 1-based seq name convention.
            // This section's a bit hairy; worked by trial and error. Final indices are respective to
            // the non-complementary seq, for both forward and reverse reading frames.
            let range = match orf {
                ReadingFrame::Fwd0 | ReadingFrame::Fwd1 | ReadingFrame::Fwd2 => {
                    (frame_open.unwrap() + 1 + offset, i + 3 + offset)
                }
                _ => (
                    seq_len_full - (i + 2 + offset),
                    seq_len_full - (frame_open.unwrap() + offset),
                ),
            };

            result.push(ReadingFrameMatch { frame: orf, range });
            frame_open = None;
        }
    }

    result
}
