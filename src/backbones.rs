//! This module  contains data structures and related for vector backbones, used for cloning. It also pulls data from Addgene.

use std::{fmt, fmt::Formatter};

use strum_macros::EnumIter;

use crate::{
    autocloning::RBS_BUFFER,
    sequence::{seq_from_str, Feature, SeqTopology},
    util::RangeIncl,
    Seq,
};

#[derive(Clone, Copy, PartialEq, EnumIter)]
pub enum AntibioticResistance {
    /// Or carbenecillin
    Ampicillin,
    Kanamycin,
    // todo: A/R
}

impl fmt::Display for AntibioticResistance {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let val = match self {
            Self::Ampicillin => "Ampicillin",
            Self::Kanamycin => "Kanamycin",
        };

        write!(f, "{val}")
    }
}

#[derive(Clone, Copy, PartialEq, EnumIter)]
pub enum ExpressionHost {
    Bacterial,
    Yeast,
    Mamallian,
}

impl fmt::Display for ExpressionHost {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let val = match self {
            Self::Bacterial => "Bacterial",
            Self::Yeast => "Yeast",
            Self::Mamallian => "Mamallian",
        };

        write!(f, "{val}")
    }
}

#[derive(Clone, Copy, PartialEq, EnumIter)]
pub enum CopyNumber {
    Unknown,
    High,
    Low,
}

impl fmt::Display for CopyNumber {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let val = match self {
            Self::Unknown => "Unknown",
            Self::High => "High",
            Self::Low => "Low",
        };

        write!(f, "{val}")
    }
}

#[allow(non_camel_case_types)]
#[derive(Clone, Copy, PartialEq, EnumIter)]
pub enum ExpressionSystem {
    None,
    T7,
    pGex,
    pBad,
    // todo A/R
}

impl fmt::Display for ExpressionSystem {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let val = match self {
            Self::None => "None",
            Self::T7 => "T7",
            Self::pGex => "pGex",
            Self::pBad => "pBad",
        };

        write!(f, "{val}")
    }
}

impl Default for ExpressionSystem {
    fn default() -> Self {
        Self::None
    }
}

#[allow(non_camel_case_types)]
#[derive(Clone, Copy, PartialEq)]
pub enum BackboneGene {
    Gfp,
    Luciferase,
    // todo A/R
}

pub struct Backbone {
    pub name: String,
    pub addgene_id: Option<u32>,
    pub seq: Seq,
    pub features: Vec<Feature>,
    pub promoter: RangeIncl, // todo: DirectioN?
    pub rbs: RangeIncl,      // todo: Features for these regions?
    pub antibiotic_resistance: AntibioticResistance,
    // pub expression_hosts: Vec<ExpressionHost>, // todo: For  now, bacterial for all
    // todo: Data on ORI and related,
    // todo:  tags etc
    pub expression_system: ExpressionSystem,
    pub genes_included: Vec<BackboneGene>,
    pub copy_number: CopyNumber,
    /// Note: This only refers to a His tag on the vector; not the insert.
    pub his_tag: Option<RangeIncl>,
    pub seq_topology: SeqTopology,
}

#[derive(Clone, Copy, PartialEq)]
pub enum CloningTechnique {
    RestrictionEnzyme,
    Pcr,
}

impl Backbone {
    /// Find the vector insertion point.
    pub fn insert_loc(&self, technique: CloningTechnique) -> Option<usize> {
        match technique {
            CloningTechnique::RestrictionEnzyme => None,
            CloningTechnique::Pcr => {
                //Initial hack: A fixed distance downstream of the RBS.
                // todo: Next: Verify in-frame.

                Some(self.rbs.end + RBS_BUFFER)
            }
        }
    }

    pub fn addgene_url(&self) -> Option<String> {
        self.addgene_id
            .map(|id| format!("https://www.addgene.org/{id}/"))
    }
}

#[derive(Default)]
pub struct BackboneFilters {
    pub host: Option<ExpressionHost>,
    pub antibiotic_resistance: Option<AntibioticResistance>,
    pub expression_system: Option<ExpressionSystem>,
    pub copy_number: Option<CopyNumber>,
    pub his_tagged: bool,
}

impl<'a> BackboneFilters {
    pub fn apply(&self, backbones: &'a [Backbone]) -> Vec<&'a Backbone> {
        let mut result = Vec::new();

        for bb in backbones {
            // todo: A/R
            if let Some(host) = self.host {}

            if let Some(es) = self.expression_system {
                if es != bb.expression_system {
                    continue;
                }
            }

            if let Some(es) = self.antibiotic_resistance {
                if es != bb.antibiotic_resistance {
                    continue;
                }
            }

            if let Some(cn) = self.copy_number {
                if cn != bb.copy_number {
                    continue;
                }
            }

            if self.his_tagged && !bb.his_tag.is_some() {
                continue;
            }

            result.push(bb);
        }

        result
    }
}

/// Load this at program init.
pub fn load_backbone_library() -> Vec<Backbone> {
    vec![
        Backbone {
            name: "pET His6 LIC cloning vector (2Bc-T)".to_owned(),
            addgene_id: Some(37_236),
            // seq: seq_from_str(""),
            seq: seq_from_str("ttcttgaagacgaaagggcctcgtgatacgcctatttttataggttaatgtcatgataataatggtttcttagacgtcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataacattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtgttgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgcagcaatggcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagccctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataatctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgtccttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcatatatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagtatacactccgctatcgctacgtgactgggtcatggctgcgccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgaggcagctgcggtaaagctcatcagcgtggtcgtgaagcgattcacagatgtctgcctgttcatccgcgtccagctcgttgagtttctccagaagcgttaatgtctggcttctgataaagcgggccatgttaagggcggttttttcctgtttggtcactgatgcctccgtgtaagggggatttctgttcatgggggtaatgataccgatgaaacgagagaggatgctcacgatacgggttactgatgatgaacatgcccggttactggaacgttgtgagggtaaacaactggcggtatggatgcggcgggaccagagaaaaatcactcagggtcaatgccagcgcttcgttaatacagatgtaggtgttccacagggtagccagcagcatcctgcgatgcagatccggaacataatggtgcagggcgctgacttccgcgtttccagactttacgaaacacggaaaccgaagaccattcatgttgttgctcaggtcgcagacgttttgcagcagcagtcgcttcacgttcgctcgcgtatcggtgattcattctgctaaccagtaaggcaaccccgccagcctagccgggtcctcaacgacaggagcacgatcatgcgcacccgtggccaggacccaacgctgcccgagatgcgccgcgtgcggctgctggagatggcggacgcgatggatatgttctgccaagggttggtttgcgcattcacagttctccgcaagaattgattggctccaattcttggagtggtgaatccgttagcgaggtgccgccggcttccattcaggtcgaggtggcccggctccatgcaccgcgacgcaacgcggggaggcagacaaggtatagggcggcgcctacaatccatgccaacccgttccatgtgctcgccgaggcggcataaatcgccgtgacgatcagcggtccagtgatcgaagttaggctggtaagagccgcgagcgatccttgaagctgtccctgatggtcgtcatctacctgcctggacagcatggcctgcaacgcgggcatcccgatgccgccggaagcgagaagaatcataatggggaaggccatccagcctcgcgtcgcgaacgccagcaagacgtagcccagcgcgtcggccgccatgccggcgataatggcctgcttctcgccgaaacgtttggtggcgggaccagtgacgaaggcttgagcgagggcgtgcaagattccgaataccgcaagcgacaggccgatcatcgtcgcgctccagcgaaagcggtcctcgccgaaaatgacccagagcgctgccggcacctgtcctacgagttgcatgataaagaagacagtcataagtgcggcgacgatagtcatgccccgcgcccaccggaaggagctgactgggttgaaggctctcaagggcatcggtcgacgctctcccttatgcgactcctgcattaggaagcagcccagtagtaggttgaggccgttgagcaccgccgccgcaaggaatggtgcatgcaaggagatggcgcccaacagtcccccggccacggggcctgccaccatacccacgccgaaacaagcgctcatgagcccgaagtggcgagcccgatcttccccatcggtgatgtcggcgatataggcgccagcaaccgcacctgtggcgccggtgatgccggccacgatgcgtccggcgtagaggatcgagatctcgatcccgcgaaattaatacgactcactatagggagaccacaacggtttccctctagtgccggctccggagagctctttaattaagcggccgccctgcaggactcgagttctagaaataattttgtttaactttaagaaggagatatagttaacctctacttccaatccggctctcatcaccatcaccatcactaataaccaactccataaggatccgcgatcgcggcgcgccacctggtggccggccggtaccacgcgtgcgcgctgatccggctgctaacaaagcccgaaaggaagctgagttggctgctgccaccgctgagcaataactagcataaccccttggggcctctaaacgggtcttgaggggttttttgctgaaaggaggaactatatccggacatccacaggacgggtgtggtcgccatgatcgcgtagtcgatagtggctccaagtagcgaagcgagcaggactgggcggcggccaaagcggtcggacagtgctccgagaacgggtgcgcatagaaattgcatcaacgcatatagcgctagcagcacgccatagtgactggcgatgctgtcggaatggacgacatcccgcaagaggcccggcagtaccggcataaccaagcctatgcctacagcatccagggtgacggtgccgaggatgacgatgagcgcattgttagatttcatacacggtgcctgactgcgttagcaatttaactgtgataaactaccgcattaaagcttatcgatgataagctgtcaaacatgagaa"),
            features: Vec::new(),
            promoter: RangeIncl::new(4010, 4029),
            rbs: RangeIncl::new(4117, 4139),
            antibiotic_resistance: AntibioticResistance::Ampicillin,
            expression_system: ExpressionSystem::T7,
            genes_included: Vec::new(),
            copy_number: CopyNumber::Low,
            his_tag: Some(RangeIncl::new(4171, 4188)),
            seq_topology: SeqTopology::Circular,
        },
        Backbone {
            name: "pET LIC cloning vector (2A-T)".to_owned(),
            addgene_id: Some(29_665),
            seq: seq_from_str("ttcttgaagacgaaagggcctcgtgatacgcctatttttataggttaatgtcatgataataatggtttcttagacgtcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataacattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtgttgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgcagcaatggcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagccctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataatctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgtccttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcatatatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagtatacactccgctatcgctacgtgactgggtcatggctgcgccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgaggcagctgcggtaaagctcatcagcgtggtcgtgaagcgattcacagatgtctgcctgttcatccgcgtccagctcgttgagtttctccagaagcgttaatgtctggcttctgataaagcgggccatgttaagggcggttttttcctgtttggtcactgatgcctccgtgtaagggggatttctgttcatgggggtaatgataccgatgaaacgagagaggatgctcacgatacgggttactgatgatgaacatgcccggttactggaacgttgtgagggtaaacaactggcggtatggatgcggcgggaccagagaaaaatcactcagggtcaatgccagcgcttcgttaatacagatgtaggtgttccacagggtagccagcagcatcctgcgatgcagatccggaacataatggtgcagggcgctgacttccgcgtttccagactttacgaaacacggaaaccgaagaccattcatgttgttgctcaggtcgcagacgttttgcagcagcagtcgcttcacgttcgctcgcgtatcggtgattcattctgctaaccagtaaggcaaccccgccagcctagccgggtcctcaacgacaggagcacgatcatgcgcacccgtggccaggacccaacgctgcccgagatgcgccgcgtgcggctgctggagatggcggacgcgatggatatgttctgccaagggttggtttgcgcattcacagttctccgcaagaattgattggctccaattcttggagtggtgaatccgttagcgaggtgccgccggcttccattcaggtcgaggtggcccggctccatgcaccgcgacgcaacgcggggaggcagacaaggtatagggcggcgcctacaatccatgccaacccgttccatgtgctcgccgaggcggcataaatcgccgtgacgatcagcggtccagtgatcgaagttaggctggtaagagccgcgagcgatccttgaagctgtccctgatggtcgtcatctacctgcctggacagcatggcctgcaacgcgggcatcccgatgccgccggaagcgagaagaatcataatggggaaggccatccagcctcgcgtcgcgaacgccagcaagacgtagcccagcgcgtcggccgccatgccggcgataatggcctgcttctcgccgaaacgtttggtggcgggaccagtgacgaaggcttgagcgagggcgtgcaagattccgaataccgcaagcgacaggccgatcatcgtcgcgctccagcgaaagcggtcctcgccgaaaatgacccagagcgctgccggcacctgtcctacgagttgcatgataaagaagacagtcataagtgcggcgacgatagtcatgccccgcgcccaccggaaggagctgactgggttgaaggctctcaagggcatcggtcgacgctctcccttatgcgactcctgcattaggaagcagcccagtagtaggttgaggccgttgagcaccgccgccgcaaggaatggtgcatgcaaggagatggcgcccaacagtcccccggccacggggcctgccaccatacccacgccgaaacaagcgctcatgagcccgaagtggcgagcccgatcttccccatcggtgatgtcggcgatataggcgccagcaaccgcacctgtggcgccggtgatgccggccacgatgcgtccggcgtagaggatcgagatctcgatcccgcgaaattaatacgactcactatagggagaccacaacggtttccctctagtgccggctccggagagctctttaattaagcggccgccctgcaggactcgagttctagaaataattttgtttaactttaagaaggagatatagatatcccaactccataaggatccgcgatcgcggcgcgccacctggtggccggccggtaccacgcgtgcgcgctgatccggctgctaacaaagcccgaaaggaagctgagttggctgctgccaccgctgagcaataactagcataaccccttggggcctctaaacgggtcttgaggggttttttgctgaaaggaggaactatatccggacatccacaggacgggtgtggtcgccatgatcgcgtagtcgatagtggctccaagtagcgaagcgagcaggactgggcggcggccaaagcggtcggacagtgctccgagaacgggtgcgcatagaaattgcatcaacgcatatagcgctagcagcacgccatagtgactggcgatgctgtcggaatggacgacatcccgcaagaggcccggcagtaccggcataaccaagcctatgcctacagcatccagggtgacggtgccgaggatgacgatgagcgcattgttagatttcatacacggtgcctgactgcgttagcaatttaactgtgataaactaccgcattaaagcttatcgatgataagctgtcaaacatgagaa"),
            features: Vec::new(), // todo temp
            promoter: RangeIncl::new(4010, 4029),
            rbs: RangeIncl::new(417, 4139),
            antibiotic_resistance: AntibioticResistance::Ampicillin,
            expression_system: ExpressionSystem::T7,
            genes_included: Vec::new(),
            copy_number: CopyNumber::Low,
            his_tag: None,
            seq_topology: SeqTopology::Circular,
        },
    ]
}
