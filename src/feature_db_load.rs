//! A library of known sequences we can use to automatically add features to a sequence.
//! This will eventually load data from an online source. Currently, it includes some common items,
//! and is stored in memory.
//!
//! Some common genes: https://en.vectorbuilder.com/search/gene.html

use crate::{
    amino_acids::AminoAcid,
    sequence::{
        seq_complement, seq_from_str, Feature,
        FeatureType::{
            self, AntibioticResistance, CodingRegion, Ori, Promoter, ProteinBind, RibosomeBindSite,
            Terminator,
        },
        Nucleotide,
        Nucleotide::{A, G, T},
        ReadingFrame, ReadingFrameMatch, Seq,
    },
    util::{match_subseq, RangeIncl},
};

struct FeatureMapItem {
    name: String,
    feature_type: FeatureType,
    seq: Seq,
}

impl FeatureMapItem {
    pub fn new(name: &str, feature_type: FeatureType, seq: Seq) -> Self {
        Self {
            name: name.to_string(),
            feature_type,
            seq,
        }
    }
}

/// Find 6x+ HIS tags in a sequence
/// todo: Currently very similar to find_reading_frame_matches
pub fn find_his_tags(seq: &[Nucleotide]) -> Vec<Feature> {
    let mut result = Vec::new();

    let seq_len_full = seq.len();
    if seq_len_full < 3 {
        return result;
    }

    for orf in [
        ReadingFrame::Fwd0,
        ReadingFrame::Fwd1,
        ReadingFrame::Fwd2,
        ReadingFrame::Rev0,
        ReadingFrame::Rev1,
        ReadingFrame::Rev2,
    ] {
        let mut offset = orf.offset();

        let seq_ = &match orf {
            ReadingFrame::Fwd0 | ReadingFrame::Fwd1 | ReadingFrame::Fwd2 => seq.to_vec(),
            _ => seq_complement(seq),
        }[offset..];

        let len = seq_.len();

        let mut tag_open = None; // Inner: Start index.

        for i_ in 0..len / 3 {
            let i = i_ * 3; // The actual sequence index.

            let nts = &seq_[i..i + 3];

            let mut matched = false;
            if let Some(aa) = AminoAcid::from_codons(nts.try_into().unwrap()) {
                if aa == AminoAcid::His {
                    matched = true;
                }
            }

            if tag_open.is_none() && matched {
                tag_open = Some(i);
            } else if tag_open.is_some() && !matched {
                let his_len = (i - tag_open.unwrap()) / 3;

                if his_len >= 6 {
                    let range = match orf {
                        ReadingFrame::Fwd0 | ReadingFrame::Fwd1 | ReadingFrame::Fwd2 => {
                            RangeIncl::new(tag_open.unwrap() + 1 + offset, i + 2 + offset)
                        }
                        _ => RangeIncl::new(
                            seq_len_full - (i + 2 + offset),
                            seq_len_full - (tag_open.unwrap() + offset) - 1,
                        ),
                    };

                    result.push(Feature {
                        range,
                        feature_type: CodingRegion,
                        label: format!("{his_len}×His"),
                        ..Default::default()
                    });
                }
                tag_open = None;
            }
        }
    }
    result
}

/// Find common promoters and Oris.
fn find_promoters(seq: &[Nucleotide]) -> Vec<Feature> {
    let items = vec![
        FeatureMapItem::new("AmpR promoter", Promoter, seq_from_str(
            "cgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataacattgaaaaaggaagagt")),

        // Perhaps slightly different from the above AmpR promoter?
        FeatureMapItem::new("AmpR promoter", Promoter, seq_from_str(
            "cgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagt")),


        FeatureMapItem::new("T3 promoter", Promoter, seq_from_str(
            "aattaaccctcactaaagg")),

        FeatureMapItem::new("lac promoter", Promoter, seq_from_str("tttacactttatgcttccggctcgtatgttg")),
        FeatureMapItem::new("lac operator", ProteinBind, seq_from_str("ttgtgagcggataacaa")),

        FeatureMapItem::new("T7 promoter", Promoter, seq_from_str("taatacgactcactatagg")),
        FeatureMapItem::new("RBS", RibosomeBindSite, seq_from_str("tttgtttaactttaagaaggaga")),
        FeatureMapItem::new("T7 terminator", Terminator, seq_from_str("ctagcataaccccttggggcctctaaacgggtcttgaggggttttttg")),

        FeatureMapItem::new("SV40 promoter", Promoter, seq_from_str("ctgaggcggaaagaaccagctgtggaatgtgtgtcagttagggtgtggaaagtccccaggctccccagcaggcagaagtatgcaaagcatgcatctcaattagtcagcaaccaggtgtggaaagtccccaggctccccagcaggcagaagtatgcaaagcatgcatctcaattagtcagcaaccatagtcccgcccctaactccgcccatcccgcccctaactccgcccagttccgcccattctccgccccatggctgactaattttttttatttatgcagaggccgaggccgcctcggcctctgagctattccagaagtagtgaggaggcttttttggaggcctaggcttttgcaaa")),
        FeatureMapItem::new("SV40 ori", Ori, seq_from_str("atcccgcccctaactccgcccagttccgcccattctccgccccatggctgactaattttttttatttatgcagaggccgaggccgcctcggcctctgagctattccagaagtagtgaggaggcttttttggaggcc")),

        FeatureMapItem::new("P15A ori", Ori, seq_from_str("ttttccataggctccgcccccctgacaagcatcacgaaatctgacgctcaaatcagtggtggcgaaacccgacaggactataaagataccaggcgtttccccctggcggctccctcgtgcgctctcctgttcctgcctttcggtttaccggtgtcattccgctgttatggccgcgtttgtctcattccacgcctgacactcagttccgggtaggcagttcgctccaagctggactgtatgcacgaaccccccgttcagtccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggaaagacatgcaaaagcaccactggcagcagccactggtaattgatttagaggagttagtcttgaagtcatgcgccggttaaggctaaactgaaaggacaagttttggtgactgcgctcctccaagccagttacctcggttcaaagagttggtagctcagagaaccttcgaaaaaccgccctgcaaggcggttttttcgttttcagagcaagagattacgcgcagaccaaaacgatctcaa")),

        FeatureMapItem::new("CAP binding site", ProteinBind, seq_from_str("atgagtgagctaactcacatta")),
        FeatureMapItem::new("pHluorin2", CodingRegion, seq_from_str("ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGAGCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACGAGCACCTGGTGTACATCATGGCCGACAAGCAGAAGAACGGCACCAAGGCCATCTTCCAGGTGCACCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGCACACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCACGGCATGGACGAGCTGTACAAGTAA")),

        // For now, we also include longer common seqs here.
        FeatureMapItem::new("f1 ori", Ori, seq_from_str("acgcgccctgtagcggcgcattaagcgcggcgggtgtggtggttacgcgcagcgtgaccgctacacttgccagcgccctagcgcccgctcctttcgctttcttcccttcctttctcgccacgttcgccggctttccccgtcaagctctaaatcgggggctccctttagggttccgatttagtgctttacggcacctcgaccccaaaaaacttgattagggtgatggttcacgtagtgggccatcgccctgatagacggtttttcgccctttgacgttggagtccacgttctttaatagtggactcttgttccaaactggaacaacactcaaccctatctcggtctattcttttgatttataagggattttgccgatttcggcctattggttaaaaaatgagctgatttaacaaaaatttaacgcgaattttaacaaaatattaacgcttacaattt")),

        FeatureMapItem::new("AmpR", AntibioticResistance, seq_from_str("gttaccaatgcttaatcagtgaggcacctatctcagcgatctgtctatttcgttcatccatagttgcctgactccccgtcgtgtagataactacgatacgggagggcttaccatctggccccagtgctgcaatgataccgcgagacccacgctcaccggctccagatttatcagcaataaaccagccagccggaagggccgagcgcagaagtggtcctgcaactttatccgcctccatccagtctattaattgttgccgggaagctagagtaagtagttcgccagttaatagtttgcgcaacgttgttgccattgctacaggcatcgtggtgtcacgctcgtcgtttggtatggcttcattcagctccggttcccaacgatcaaggcgagttacatgatcccccatgttgtgcaaaaaagcggttagctccttcggtcctccgatcgttgtcagaagtaagttggccgcagtgttatcactcatggttatggcagcactgcataattctcttactgtcatgccatccgtaagatgcttttctgtgactggtgagtactcaaccaagtcattctgagaatagtgtatgcggcgaccgagttgctcttgcccggcgtcaatacgggataataccgcgccacatagcagaactttaaaagtgctcatcattggaaaacgttcttcggggcgaaaactctcaaggatcttaccgctgttgagatccagttcgatgtaacccactcgtgcacccaactgatcttcagcatcttttactttcaccagcgtttctgggtgagcaaaaacaggaaggcaaaatgccgcaaaaaagggaataagggcgacacggaaatgttgaatactcat")),
        FeatureMapItem::new("AmpR", AntibioticResistance, seq_from_str("atgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgtagcaatggcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtcgcgcggtatcattgcagcactggggccagatggtaagccctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaa")),
        FeatureMapItem::new("KanR", AntibioticResistance, seq_from_str("gttagaaaaactcatcgagcatcaaatgaaactgcaatttattcatatcaggattatcaataccatatttttgaaaaagccgtttctgtaatgaaggagaaaactcaccgaggcagttccataggatggcaagatcctggtatcggtctgcgattccgactcgtccaacatcaatacaacctattaatttcccctcgtcaaaaataaggttatcaagtgagaaatcaccatgagtgacgactgaatccggtgagaatggcaaaagtttatgcatttctttccagacttgttcaacaggccagccattacgctcgtcatcaaaatcactcgcatcaaccaaaccgttattcattcgtgattgcgcctgagcgagacgaaatacgcgatcgctgttaaaaggacaattacaaacaggaatcgaatgcaaccggcgcaggaacactgccagcgcatcaacaatattttcacctgaatcaggatattcttctaatacctggaatgctgttttcccagggatcgcagtggtgagtaaccatgcatcatcaggagtacggataaaatgcttgatggtcggaagaggcataaattccgtcagccagtttagtctgaccatctcatctgtaacatcattggcaacgctacctttgccatgtttcagaaacaactctggcgcatcgggcttcccatacaatcgatagattgtcgcacctgattgcccgacattatcgcgagcccatttatacccatataaatcagcatccatgttggaatttaatcgcggcctagagcaagacgtttcccgttgaatatggctcat")),
        FeatureMapItem::new("CmR", AntibioticResistance, seq_from_str("attacgccccgccctgccactcatcgcagtactgttgtaattcattaagcattctgccgacatggaagccatcacaaacggcatgatgaacctgaatcgccagcggcatcagcaccttgtcgccttgcgtataatatttgcccatagtgaaaacgggggcgaagaagttgtccatattggccacgtttaaatcaaaactggtgaaactcacccagggattggctgagacgaaaaacatattctcaataaaccctttagggaaataggccaggttttcaccgtaacacgccacatcttgcgaatatatgtgtagaaactgccggaaatcgtcgtggtattcactccagagcgatgaaaacgtttcagtttgctcatggaaaacggtgtaacaagggtgaacactatcccatatcaccagctcaccgtctttcattgccatacggaactccggatgagcattcatcaggcgggcaagaatgtgaataaaggccggataaaacttgtgcttatttttctttacggtctttaaaaaggccgtaatatccagctgaacggtctggttataggtacattgagcaactgactgaaatgcctcaaaatgttctttacgatgccattgggatatatcaacggtggtatatccagtgatttttttctccat")),

        FeatureMapItem::new("EGFP", CodingRegion, seq_from_str("ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA")),

        // todo: How general is this?
        FeatureMapItem::new("Ori", Ori, seq_from_str("ttttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaagaacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaa")),

    ];

    let mut result = Vec::new();

    for item in items {
        let (matches_fwd, matches_rev) = match_subseq(&item.seq, seq);

        // todo: Consider short descriptive notes for each.
        for range in matches_fwd {
            result.push(Feature {
                range,
                feature_type: item.feature_type,
                label: item.name.clone(),
                // todo: Direction?
                ..Default::default()
            });
        }

        let seq_len = seq.len();
        for range in matches_rev {
            result.push(Feature {
                // todo: QC off-by-one errors etc.
                range: RangeIncl::new(seq_len - range.end, seq_len - range.start),
                feature_type: item.feature_type,
                label: item.name.clone(),
                // todo: Direction?
                ..Default::default()
            });
        }
    }

    result
}

pub fn find_features(seq: &[Nucleotide]) -> Vec<Feature> {
    let mut result = Vec::new();

    result.append(&mut find_his_tags(seq));
    result.append(&mut find_promoters(seq));

    result
}
