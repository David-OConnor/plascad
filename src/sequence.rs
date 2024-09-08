use std::{fmt::Display, io};

use bincode::{Decode, Encode};
use chrono::NaiveDate;
use num_enum::TryFromPrimitive;

use crate::{
    primer::PrimerDirection,
    sequence::Nucleotide::{A, C, G, T},
    util::{match_subseq, RangeIncl},
    Color,
};
pub const MIN_SEARCH_LEN: usize = 3;

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
    T = 0b00,
    C = 0b01,
    A = 0b10,
    G = 0b11,
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

    pub fn complement(self) -> Self {
        match self {
            A => T,
            T => A,
            G => C,
            C => G,
        }
    }

    /// Molecular weight, in Daltons, in a DNA strand.
    /// [Weight source: NorthWestern](http://biotools.nubic.northwestern.edu/OligoCalc.html)
    pub fn weight(&self) -> f32 {
        match self {
            A => 313.21,
            T => 304.2,
            G => 329.21,
            C => 289.18,
        }
    }

    /// Optical density of a 1mL solution, in a cuvette with 1cm pathlength.
    /// Result is in nm.
    /// http://biotools.nubic.northwestern.edu/OligoCalc.html
    pub fn a_max(&self) -> f32 {
        match self {
            A => 259.,
            T => 267.,
            G => 253.,
            C => 271.,
        }
    }

    /// Optical density of a 1mL solution, in a cuvette with 1cm pathlength.
    /// Result is in 1/(Moles x cm)
    /// http://biotools.nubic.northwestern.edu/OligoCalc.html
    pub fn molar_density(&self) -> f32 {
        match self {
            A => 15_200.,
            T => 8_400.,
            G => 12_010.,
            C => 7_050.,
        }
    }
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
    /// Like Primer, this is not a real feature; we use it to draw the selection highlighted area.
    Selection,
    /// Ie operators.
    ProteinBind,
    Terminator,
}

impl Default for FeatureType {
    fn default() -> Self {
        Self::Generic
    }
}

impl FeatureType {
    /// For displaying in the UI
    pub fn to_string(self) -> String {
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
            Self::Selection => "",
            Self::ProteinBind => "Operator",
            Self::Terminator => "Terminator",
        }
        .to_owned()
    }

    pub fn color(&self) -> Color {
        match self {
            Self::Generic => (255, 0, 255),
            Self::Gene => (255, 128, 128),
            Self::Ori => (40, 200, 128),
            Self::RibosomeBindSite => (255, 204, 252),
            Self::Promoter => (240, 190, 70),
            Self::AntibioticResistance => (0, 200, 110),
            Self::Primer => (0, 0, 0),             // N/A for now at least.
            Self::CodingRegion => (100, 200, 255), // N/A for now at least.
            Self::LongTerminalRepeat => (150, 200, 255), // N/A for now at least.
            Self::Source => (120, 70, 120),
            Self::Exon => (255, 255, 180),
            Self::Transcript => (180, 255, 180),
            Self::Selection => (255, 255, 0),
            Self::ProteinBind => (128, 110, 150),
            Self::Terminator => (255, 110, 150),
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
            "protein_bind" => Self::ProteinBind,
            "terminator" => Self::Terminator,
            _ => Self::Generic,
        }
    }

    /// Create a string for use with SnapGene and GenBank formats.
    pub fn to_external_str(self) -> String {
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
            Self::Selection => "",
            Self::ProteinBind => "protein_bind",
            Self::Terminator => "terminator",
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
    pub fn to_string(self) -> String {
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
    // pub range: (usize, usize),
    /// 1-based indexing, inclusive. (Note: Could also use the builtin RangeInclusive.)
    pub range: RangeIncl,
    pub feature_type: FeatureType,
    pub direction: FeatureDirection,
    pub label: String,
    /// By default, we display features using featuretype-specific color. Allow the user
    /// to override this.
    pub color_override: Option<Color>,
    pub notes: Vec<(String, String)>,
}

impl Default for Feature {
    fn default() -> Self {
        Self {
            range: RangeIncl::new(1, 1),
            feature_type: Default::default(),
            direction: Default::default(),
            label: Default::default(),
            color_override: Default::default(),
            notes: Default::default(),
        }
    }
}

impl Feature {
    pub fn label(&self) -> String {
        if self.label.is_empty() {
            self.feature_type.to_string()
        } else {
            self.label.clone()
        }
    }

    /// Get the color to draw; type color, unless overridden.
    pub fn color(&self) -> Color {
        match self.color_override {
            Some(c) => c,
            None => self.feature_type.color(),
        }
    }

    /// Get the feature len, in usize.
    pub fn len(&self, seq_len: usize) -> usize {
        if self.range.end > self.range.start {
            self.range.end - self.range.start + 1
        } else {
            // ie a wrap through the origin
            self.range.end + seq_len - self.range.start + 1
        }
    }

    /// Formats the indexes, and size of this feature.
    pub fn location_descrip(&self, seq_len: usize) -> String {
        format!(
            "{}..{}  {} bp",
            self.range.start,
            self.range.end,
            self.len(seq_len)
        )
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
        *nt = nt.complement();
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

/// Contains sequence-level metadata.
#[derive(Clone, Default, Encode, Decode)]
pub struct Metadata {
    pub plasmid_name: String,
    pub comments: Vec<String>,
    pub references: Vec<Reference>,
    pub locus: String,
    // pub date: Option<NaiveDate>,
    /// We use this to prevent manual encode/decode for the whoel struct.
    pub date: Option<(i32, u8, u8)>, // Y M D
    pub definition: Option<String>,
    pub molecule_type: Option<String>,
    pub division: String,
    pub accession: Option<String>,
    pub version: Option<String>,
    // pub keywords: Vec<String>,
    pub keywords: Option<String>, // todo vec?
    pub source: Option<String>,
    pub organism: Option<String>,
}

/// Based on GenBank's reference format
#[derive(Default, Clone, Encode, Decode)]
pub struct Reference {
    pub description: String,
    pub authors: Option<String>,
    pub consortium: Option<String>,
    pub title: String,
    pub journal: Option<String>,
    pub pubmed: Option<String>,
    pub remark: Option<String>,
}

pub struct SearchMatch {
    /// 0-based indexing.
    pub range: RangeIncl,
    // todo: More A/R
}

/// Find exact matches in the target sequence of our search nucleotides.
/// todo: Optionally support partial matches.
/// todo:
pub fn find_search_matches(seq: &[Nucleotide], search_seq: &[Nucleotide]) -> Vec<SearchMatch> {
    let (mut fwd, mut rev) = match_subseq(search_seq, seq);

    fwd.append(&mut rev);
    fwd.into_iter().map(|range| SearchMatch { range }).collect()
}

/// Sequence weight, in Daltons. Assumes single-stranded.
pub fn seq_weight(seq: &[Nucleotide]) -> f32 {
    let mut result = 0.;

    for nt in seq {
        result += nt.weight();
    }

    // result -= WATER_WEIGHT * (seq.len() - 1) as f32;

    result -= 61.96;

    result
}
