use std::fmt;

use bincode::{Decode, Encode};
use na_seq::{Nucleotide, Nucleotide::*};

#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum AaIdent {
    OneLetter,
    ThreeLetters,
}

/// This struct and its methods are largely copied from the `peptide` project.
#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum AminoAcid {
    Arg,
    His,
    Lys,
    Asp,
    Glu,
    Ser,
    Thr,
    Asn,
    Gln,
    Cys,
    Sec,
    Gly,
    Pro,
    Ala,
    Val,
    Ile,
    Leu,
    Met,
    Phe,
    Tyr,
    Trp,
}

impl AminoAcid {
    pub fn ident_single_letter(&self) -> String {
        match self {
            Self::Arg => "R",
            Self::His => "H",
            Self::Lys => "K",
            Self::Asp => "D",
            Self::Glu => "E",
            Self::Ser => "S",
            Self::Thr => "T",
            Self::Asn => "N",
            Self::Gln => "Q",
            Self::Cys => "C",
            Self::Sec => "U",
            Self::Gly => "G",
            Self::Pro => "P",
            Self::Ala => "A",
            Self::Val => "V",
            Self::Ile => "I",
            Self::Leu => "L",
            Self::Met => "M",
            Self::Phe => "F",
            Self::Tyr => "Y",
            Self::Trp => "W",
        }
        .to_owned()
    }

    pub fn ident_3_letter(&self) -> String {
        match self {
            Self::Arg => "Arg",
            Self::His => "His",
            Self::Lys => "Lys",
            Self::Asp => "Asp",
            Self::Glu => "Glu",
            Self::Ser => "Ser",
            Self::Thr => "Thr",
            Self::Asn => "Asn",
            Self::Gln => "Gln",
            Self::Cys => "Cys",
            Self::Sec => "Sec",
            Self::Gly => "Gly",
            Self::Pro => "Pro",
            Self::Ala => "Ala",
            Self::Val => "Val",
            Self::Ile => "Ile",
            Self::Leu => "Leu",
            Self::Met => "Met",
            Self::Phe => "Phe",
            Self::Tyr => "Tyr",
            Self::Trp => "Trp",
        }
        .to_owned()
    }

    /// Return the molecular weight, in Da.
    /// Source: https://www.promega.com/resources/tools/amino-acid-chart-amino-acid-structure/
    /// todo: This table is not very precise; consider updating with a better source.
    pub fn weight(&self) -> f32 {
        match self {
            Self::Arg => 174.,
            Self::His => 155.,
            Self::Lys => 146.,
            Self::Asp => 133.,
            Self::Glu => 147.,
            Self::Ser => 105.,
            Self::Thr => 119.,
            Self::Asn => 132.,
            Self::Gln => 146.,
            Self::Cys => 121.,
            Self::Sec => 168.06,
            Self::Gly => 75.,
            Self::Pro => 115.,
            Self::Ala => 89.,
            Self::Val => 117.,
            Self::Ile => 131.,
            Self::Leu => 131.,
            Self::Met => 149.,
            Self::Phe => 165.,
            Self::Tyr => 181.,
            Self::Trp => 204.,
        }
    }

    /// Used for determining protein hydropathy. High (eg positive) values intdicate hydrophilic
    /// AAs. (Seems to not be completely true from some example checks? Some traditionally hydrophilic
    /// proteins like Proline (-1.6) and Glycine (-4) are on the list, but the very negative values
    /// are not associated with traditionally hydrophillic AAs.
    /// [Kyte, Doolittle](https://web.expasy.org/protscale/pscale/Hydropath.Doolittle.html)
    pub fn hydropathicity(&self) -> f32 {
        match self {
            Self::Arg => -4.5,
            Self::His => -3.2,
            Self::Lys => -3.9,
            Self::Asp => -3.5,
            Self::Glu => -3.5,
            Self::Ser => -0.8,
            Self::Thr => -0.7,
            Self::Asn => -3.5,
            Self::Gln => -3.5,
            Self::Cys => 2.5,
            Self::Sec => 0., // todo?
            Self::Gly => -0.4,
            Self::Pro => -1.6,
            Self::Ala => 1.8,
            Self::Val => 4.2,
            Self::Ile => 4.5,
            Self::Leu => 3.8,
            Self::Met => 1.9,
            Self::Phe => 2.8,
            Self::Tyr => -1.3,
            Self::Trp => -0.9,
        }
    }

    /// https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables#/media/File:Aminoacids_table.svg
    /// If a codon has less than 3 nucleotides, it means the third can be any; this may have both conciseness,
    /// and performance advantages.
    pub fn codons(&self) -> Vec<Vec<Nucleotide>> {
        match self {
            // todo: Should we do wildcards etc, to speed up matching? Ie Arg is just [C, G].
            Self::Arg => vec![vec![C, G]],
            Self::Gln => vec![vec![C, A, G], vec![C, A, A]],
            Self::His => vec![vec![C, A, C], vec![C, A, T]],
            Self::Pro => vec![vec![C, C]],
            Self::Leu => vec![vec![C, T]],
            Self::Met => vec![vec![A, T, G]],
            _ => Vec::new(),
        }
    }

    pub fn from_codons(codons: [Nucleotide; 3]) -> Option<Self> {
        // Handle cases that are defined entirely by the first two codons.
        match codons[0..2] {
            [C, G] => return Some(Self::Arg),
            [C, C] => return Some(Self::Pro),
            [C, T] => return Some(Self::Leu),
            [T, C] => return Some(Self::Ser),
            [G, G] => return Some(Self::Gly),
            [G, C] => return Some(Self::Ala),
            [G, T] => return Some(Self::Val),
            [A, C] => return Some(Self::Thr),
            _ => (),
        }

        match codons {
            [A, T, G] => Some(Self::Met),
            [A, T, A] => Some(Self::Ile),
            [A, T, C] => Some(Self::Ile),
            [A, T, T] => Some(Self::Ile),
            [C, A, G] => Some(Self::Gln),
            [C, A, A] => Some(Self::Gln),
            [C, A, C] => Some(Self::His),
            [C, A, T] => Some(Self::His),
            [T, G, G] => Some(Self::Trp),
            [T, G, A] => None,
            [T, G, C] => Some(Self::Cys),
            [T, G, T] => Some(Self::Cys),
            [T, A, G] => None,
            [T, A, A] => None,
            [T, A, C] => Some(Self::Tyr),
            [T, A, T] => Some(Self::Tyr),
            [T, T, G] => Some(Self::Leu),
            [T, T, A] => Some(Self::Leu),
            [T, T, C] => Some(Self::Phe),
            [T, T, T] => Some(Self::Phe),
            [G, A, G] => Some(Self::Glu),
            [G, A, A] => Some(Self::Glu),
            [G, A, C] => Some(Self::Asp),
            [G, A, T] => Some(Self::Asp),
            [A, G, G] => Some(Self::Arg),
            [A, G, A] => Some(Self::Arg),
            [A, G, C] => Some(Self::Ser),
            [A, G, T] => Some(Self::Ser),
            [A, A, G] => Some(Self::Lys),
            [A, A, A] => Some(Self::Lys),
            [A, A, C] => Some(Self::Asn),
            [A, A, T] => Some(Self::Asn),
            _ => unreachable!(), // This the 2-nt pattners we handled above.
        }
    }
}

impl fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let v = match self {
            Self::Arg => "Arg (R)",
            Self::His => "His (H)",
            Self::Lys => "Lys (K)",
            Self::Asp => "Asp (D)",
            Self::Glu => "Glu (E)",
            Self::Ser => "Ser (S)",
            Self::Thr => "Thr (T)",
            Self::Asn => "Asn (N)",
            Self::Gln => "Gln (Q)",
            Self::Cys => "Cys (C)",
            Self::Sec => "Sec (U)",
            Self::Gly => "Gly (G)",
            Self::Pro => "Pro (P)",
            Self::Ala => "Ala (A)",
            Self::Val => "Val (V)",
            Self::Ile => "Ile (I)",
            Self::Leu => "Leu (L)",
            Self::Met => "Met (M)",
            Self::Phe => "Phe (F)",
            Self::Tyr => "Tyr (Y)",
            Self::Trp => "Trp (W)",
        };

        write!(f, "{}", v)
    }
}
