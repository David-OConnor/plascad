//! This module contains info related to Restriction Enzyme sites.
//!
//! [Wikipedia: List of RE sites](https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites:_A)
//! [NEB guide](https://www.neb.com/en-us/tools-and-resources/selection-charts/frequencies-of-restriction-sites)
//!
//! Note: This module only currently includes a selection of popular REs, and only ones that match
//! exact NTs.

use std::{
    collections::HashMap,
    hash::{Hash, Hasher},
};

use crate::sequence::{
    seq_complement, seq_to_str,
    Nucleotide::{self, A, C, G, T},
    Seq,
};

/// Unlike `Nucleotide`, this includes wildcards
#[derive(Clone, Copy)]
pub enum NucleotideGeneral {
    A,
    T,
    C,
    G,
    /// Any
    N,
    /// A or T
    W,
    /// C or G
    S,
    /// Pyrimidines: C or T
    Y,
    /// Purines: A or G
    R,
    /// A or C
    M,
    /// G or T
    K,
}

#[derive(Debug)]
pub struct ReMatch {
    pub lib_index: usize,
    /// Cuts after this index, in the "forward" direction.
    pub seq_index: usize,
    // pub direction: PrimerDirection,
    /// todo: Experimenting
    /// The number of matches found for this RE.
    pub match_count: usize,
}

#[derive(Clone, Eq)]
pub struct RestrictionEnzyme {
    pub name: String,
    /// From the 5' end.
    pub seq: Seq, // todo: You may eventually need Vec<NucleotideGeneral>.
    /// Index to cut after, from the 5' end. For blunt ends, this will be
    /// halfway through the seq (rounded down)
    pub cut_after: u8,
}

impl Hash for RestrictionEnzyme {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.name.hash(state);
    }
}

impl PartialEq for RestrictionEnzyme {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
    }
}

impl RestrictionEnzyme {
    pub fn new(name: &str, seq: Seq, cut_after: u8) -> Self {
        Self {
            name: name.to_owned(),
            seq,
            cut_after,
        }
    }

    pub fn makes_blunt_ends(&self) -> bool {
        self.cut_after as isize + 1 == self.seq.len() as isize / 2
    }

    /// A depiction of where to cut.
    pub fn cut_depiction(&self) -> String {
        let mut nt_chars = seq_to_str(&self.seq);

        let mut result = String::new();

        for (i, nt_char) in nt_chars.chars().enumerate() {
            result.push(nt_char);
            if i as u8 == self.cut_after {
                result.push_str(" | ");
            }
        }

        result
    }

    /// Find the overhanging NTs 5' of a sequence's top strand.
    pub fn overhang_top_left(&self) -> Vec<Nucleotide> {
        let cut = self.cut_after as usize + 1;
        let len = self.seq.len();

        if cut as isize - 2 >= len as isize / 2 {
            Vec::new() // No overhang on this strand.
        } else {
            if len - cut > cut {
                return Vec::new();
            }

            self.seq[cut..len - cut].to_vec()
        }
    }

    pub fn overhang_top_right(&self) -> Vec<Nucleotide> {
        let cut = self.cut_after as usize + 1;
        let len = self.seq.len();

        if cut as isize - 2 < len as isize / 2 {
            Vec::new() // No overhang on this strand.
        } else {
            self.seq[len - cut..cut].to_vec()
        }
    }

    pub fn overhang_bottom_left(&self) -> Vec<Nucleotide> {
        // todo: DRY implementation of reverse without compl
        let x = self.overhang_top_right();
        let mut result = x.to_vec();

        for nt in &mut result {
            *nt = nt.complement();
        }

        result
    }

    pub fn overhang_bottom_right(&self) -> Vec<Nucleotide> {
        // todo: DRY implementation of reverse without compl
        let x = self.overhang_top_left();
        let mut result = x.to_vec();

        for nt in &mut result {
            *nt = nt.complement();
        }

        result
        // seq_complement(&self.overhang_top_left())
    }
}

/// Go through a sequence, and attempt to match each lib in our RE lib to the sequence, in both directions
pub fn find_re_matches(seq: &[Nucleotide], lib: &[RestrictionEnzyme]) -> Vec<ReMatch> {
    let mut result = Vec::new();

    let mut match_counts = HashMap::new(); // lib index, count

    // todo: Reverse?

    for (lib_index, re) in lib.iter().enumerate() {
        let seq_len = seq.len();
        for i in 0..seq_len {
            if i + re.seq.len() + 1 >= seq_len {
                continue;
            }

            if re.seq == seq[i..i + re.seq.len()] {
                result.push(ReMatch {
                    lib_index,
                    // direction: PrimerDirection::Forward,
                    seq_index: i + 1, // +1 indexing.
                    match_count: 0,   // Updated below.
                });

                if match_counts.contains_key(&lib_index) {
                    *match_counts.get_mut(&lib_index).unwrap() += 1;
                } else {
                    match_counts.insert(lib_index, 1);
                }
            }
        }
    }

    // Apply match counts.
    for re_match in &mut result {
        re_match.match_count = match_counts[&re_match.lib_index];
    }

    result
}

/// Load a set of common Restriction enzymes. Call this at program start, to load into a state field.
pub fn load_re_library() -> Vec<RestrictionEnzyme> {
    vec![
        RestrictionEnzyme::new("AanI", vec![T, T, A, T, A, A], 2),
        RestrictionEnzyme::new("AatI", vec![A, G, G, C, C, T], 2),
        RestrictionEnzyme::new("AatII", vec![G, A, C, G, T, C], 4),
        RestrictionEnzyme::new("AbsI", vec![C, C, T, C, G, A, G, G], 1),
        RestrictionEnzyme::new("Acc65I", vec![G, G, T, A, C, C], 0),
        RestrictionEnzyme::new("AflII", vec![C, T, T, A, A, G], 0),
        RestrictionEnzyme::new("AgeI", vec![A, C, C, G, G, T], 0),
        RestrictionEnzyme::new("ApaI", vec![G, G, G, C, C, C], 4),
        RestrictionEnzyme::new("AscI", vec![G, G, C, G, C, G, C, C], 1),
        RestrictionEnzyme::new("AseI", vec![A, T, T, A, A, T], 1),
        RestrictionEnzyme::new("AsiSI", vec![G, C, G, A, T, C, G, C], 4),
        RestrictionEnzyme::new("BamHI", vec![G, G, A, T, C, C], 0),
        RestrictionEnzyme::new("BcII", vec![T, G, A, T, C, A], 0),
        // RestrictionEnzyme::new("BglI", vec![], 0),
        RestrictionEnzyme::new("BglII", vec![A, G, A, T, C, T], 0),
        RestrictionEnzyme::new("BmtI", vec![G, C, T, A, G, C], 4),
        RestrictionEnzyme::new("BsgDI", vec![A, T, C, G, A, T], 1),
        RestrictionEnzyme::new("BsgEI", vec![T, C, C, G, G, A], 0),
        RestrictionEnzyme::new("BsgHI", vec![T, C, A, T, G, A], 0),
        RestrictionEnzyme::new("BspEI", vec![T, C, C, G, G, A], 0),
        RestrictionEnzyme::new("BstBI", vec![T, T, C, G, A, A], 1),
        RestrictionEnzyme::new("ClaI", vec![A, T, C, G, A, T], 1),
        RestrictionEnzyme::new("EcoRI", vec![G, A, A, T, T, C], 0),
        RestrictionEnzyme::new("EcoRV", vec![G, A, T, A, T, C], 2),
        RestrictionEnzyme::new("HindIII", vec![A, A, G, C, T, T], 0),
        RestrictionEnzyme::new("FspI", vec![T, G, C, G, C, A], 2),
        // todo: TOo common
        // RestrictionEnzyme::new("HhaI", vec![G, C, G, C], 2),
        RestrictionEnzyme::new("HpaI", vec![G, T, T, A, A, C], 2),
        RestrictionEnzyme::new("KnpI", vec![G, G, T, A, C, C], 4),
        RestrictionEnzyme::new("MauBI", vec![C, G, C, G, C, G, C, G], 1),
        RestrictionEnzyme::new("MscI", vec![T, G, G, C, C, A], 2),
        RestrictionEnzyme::new("NdeI", vec![C, A, T, A, T, G], 1),
        RestrictionEnzyme::new("NotI", vec![G, C, G, G, C, C, G, C], 1),
        RestrictionEnzyme::new("NruI", vec![T, C, G, C, G, A], 2),
        RestrictionEnzyme::new("NsiI", vec![A, T, G, C, A, T], 4),
        RestrictionEnzyme::new("PacI", vec![T, T, A, A, T, T, A, A], 4),
        RestrictionEnzyme::new("PciI", vec![A, C, A, T, G, T], 0),
        RestrictionEnzyme::new("PmeI", vec![G, T, T, T, A, A, A, C], 4),
        RestrictionEnzyme::new("PmII", vec![C, A, C, G, T, G], 2),
        RestrictionEnzyme::new("PmlI", vec![C, A, C, G, T, G], 2),
        RestrictionEnzyme::new("PsiI", vec![T, T, A, T, A, A], 2),
        RestrictionEnzyme::new("PspOMI", vec![G, G, G, C, C, C], 0),
        RestrictionEnzyme::new("PstI", vec![C, T, G, C, A, G], 4),
        RestrictionEnzyme::new("SacI", vec![G, A, G, C, T, C], 4),
        // RestrictionEnzyme::new("SapI", vec![G, C, T, C, T, T, C], 4), // todo: Unclea on the cut site
        RestrictionEnzyme::new("SalI", vec![G, T, C, G, A, C], 0),
        RestrictionEnzyme::new("ScaI", vec![A, G, T, A, C, T], 2),
        RestrictionEnzyme::new("SmaI", vec![C, C, C, G, G, G], 2),
        RestrictionEnzyme::new("SfbI", vec![C, C, T, G, C, A, G, G], 5),
        RestrictionEnzyme::new("SfoI", vec![G, G, C, G, C, C], 2),
        RestrictionEnzyme::new("SpeI", vec![A, C, T, A, G, T], 0),
        RestrictionEnzyme::new("SphI", vec![G, C, A, T, G, C], 4),
        RestrictionEnzyme::new("SrfI", vec![G, C, C, C, G, G, G, C], 3),
        RestrictionEnzyme::new("StuI", vec![A, G, G, C, C, T], 2),
        RestrictionEnzyme::new("XbaI", vec![T, C, T, A, G, A], 0),
        RestrictionEnzyme::new("XhoI", vec![C, T, C, G, A, G], 0),
        RestrictionEnzyme::new("ZraI", vec![G, A, C, G, T, C], 2),
        // RestrictionEnzyme::new("HaeIII", vec![G, G, C, C], 1), // Too many matches
    ]
}
