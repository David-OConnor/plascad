//! A library of known sequences we can use to automatically add features to a sequence.
//! This will eventually load data from an online source. Currently, it includes some common items,
//! and is stored in memory.

use crate::{
    sequence::{seq_from_str, Feature, FeatureType::{self, ProteinBind, Promoter, AntibioticResistance, Ori, Terminator}, Nucleotide, Seq},
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

/// todo: Abstract this out later, if it makes sense.
fn find_his_tags(seq: &[Nucleotide]) -> Vec<Feature> {
    let mut result = Vec::new();

    result
}

// todo: Consider longer common sequences like

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
        FeatureMapItem::new("T7 terminator", Terminator, seq_from_str("ctagcataaccccttggggcctctaaacgggtcttgaggggttttttg")),

        FeatureMapItem::new("SV40 promoter", Promoter, seq_from_str("ctgaggcggaaagaaccagctgtggaatgtgtgtcagttagggtgtggaaagtccccaggctccccagcaggcagaagtatgcaaagcatgcatctcaattagtcagcaaccaggtgtggaaagtccccaggctccccagcaggcagaagtatgcaaagcatgcatctcaattagtcagcaaccatagtcccgcccctaactccgcccatcccgcccctaactccgcccagttccgcccattctccgccccatggctgactaattttttttatttatgcagaggccgaggccgcctcggcctctgagctattccagaagtagtgaggaggcttttttggaggcctaggcttttgcaaa")),
        FeatureMapItem::new("SV40 ori", Ori, seq_from_str("atcccgcccctaactccgcccagttccgcccattctccgccccatggctgactaattttttttatttatgcagaggccgaggccgcctcggcctctgagctattccagaagtagtgaggaggcttttttggaggcc")),

        FeatureMapItem::new("P15A ori", Ori, seq_from_str("ttttccataggctccgcccccctgacaagcatcacgaaatctgacgctcaaatcagtggtggcgaaacccgacaggactataaagataccaggcgtttccccctggcggctccctcgtgcgctctcctgttcctgcctttcggtttaccggtgtcattccgctgttatggccgcgtttgtctcattccacgcctgacactcagttccgggtaggcagttcgctccaagctggactgtatgcacgaaccccccgttcagtccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggaaagacatgcaaaagcaccactggcagcagccactggtaattgatttagaggagttagtcttgaagtcatgcgccggttaaggctaaactgaaaggacaagttttggtgactgcgctcctccaagccagttacctcggttcaaagagttggtagctcagagaaccttcgaaaaaccgccctgcaaggcggttttttcgttttcagagcaagagattacgcgcagaccaaaacgatctcaa")),

        FeatureMapItem::new("CAP binding site", ProteinBind, seq_from_str("atgagtgagctaactcacatta")),

        // For now, we also include longer common seqs here.
        FeatureMapItem::new("f1 ori", Ori, seq_from_str("acgcgccctgtagcggcgcattaagcgcggcgggtgtggtggttacgcgcagcgtgaccgctacacttgccagcgccctagcgcccgctcctttcgctttcttcccttcctttctcgccacgttcgccggctttccccgtcaagctctaaatcgggggctccctttagggttccgatttagtgctttacggcacctcgaccccaaaaaacttgattagggtgatggttcacgtagtgggccatcgccctgatagacggtttttcgccctttgacgttggagtccacgttctttaatagtggactcttgttccaaactggaacaacactcaaccctatctcggtctattcttttgatttataagggattttgccgatttcggcctattggttaaaaaatgagctgatttaacaaaaatttaacgcgaattttaacaaaatattaacgcttacaattt")),

        FeatureMapItem::new("AmpR", AntibioticResistance, seq_from_str("gttaccaatgcttaatcagtgaggcacctatctcagcgatctgtctatttcgttcatccatagttgcctgactccccgtcgtgtagataactacgatacgggagggcttaccatctggccccagtgctgcaatgataccgcgagacccacgctcaccggctccagatttatcagcaataaaccagccagccggaagggccgagcgcagaagtggtcctgcaactttatccgcctccatccagtctattaattgttgccgggaagctagagtaagtagttcgccagttaatagtttgcgcaacgttgttgccattgctacaggcatcgtggtgtcacgctcgtcgtttggtatggcttcattcagctccggttcccaacgatcaaggcgagttacatgatcccccatgttgtgcaaaaaagcggttagctccttcggtcctccgatcgttgtcagaagtaagttggccgcagtgttatcactcatggttatggcagcactgcataattctcttactgtcatgccatccgtaagatgcttttctgtgactggtgagtactcaaccaagtcattctgagaatagtgtatgcggcgaccgagttgctcttgcccggcgtcaatacgggataataccgcgccacatagcagaactttaaaagtgctcatcattggaaaacgttcttcggggcgaaaactctcaaggatcttaccgctgttgagatccagttcgatgtaacccactcgtgcacccaactgatcttcagcatcttttactttcaccagcgtttctgggtgagcaaaaacaggaaggcaaaatgccgcaaaaaagggaataagggcgacacggaaatgttgaatactcat")),
        FeatureMapItem::new("KanR", AntibioticResistance, seq_from_str("gttagaaaaactcatcgagcatcaaatgaaactgcaatttattcatatcaggattatcaataccatatttttgaaaaagccgtttctgtaatgaaggagaaaactcaccgaggcagttccataggatggcaagatcctggtatcggtctgcgattccgactcgtccaacatcaatacaacctattaatttcccctcgtcaaaaataaggttatcaagtgagaaatcaccatgagtgacgactgaatccggtgagaatggcaaaagtttatgcatttctttccagacttgttcaacaggccagccattacgctcgtcatcaaaatcactcgcatcaaccaaaccgttattcattcgtgattgcgcctgagcgagacgaaatacgcgatcgctgttaaaaggacaattacaaacaggaatcgaatgcaaccggcgcaggaacactgccagcgcatcaacaatattttcacctgaatcaggatattcttctaatacctggaatgctgttttcccagggatcgcagtggtgagtaaccatgcatcatcaggagtacggataaaatgcttgatggtcggaagaggcataaattccgtcagccagtttagtctgaccatctcatctgtaacatcattggcaacgctacctttgccatgtttcagaaacaactctggcgcatcgggcttcccatacaatcgatagattgtcgcacctgattgcccgacattatcgcgagcccatttatacccatataaatcagcatccatgttggaatttaatcgcggcctagagcaagacgtttcccgttgaatatggctcat")),
        FeatureMapItem::new("CmR", AntibioticResistance, seq_from_str("attacgccccgccctgccactcatcgcagtactgttgtaattcattaagcattctgccgacatggaagccatcacaaacggcatgatgaacctgaatcgccagcggcatcagcaccttgtcgccttgcgtataatatttgcccatagtgaaaacgggggcgaagaagttgtccatattggccacgtttaaatcaaaactggtgaaactcacccagggattggctgagacgaaaaacatattctcaataaaccctttagggaaataggccaggttttcaccgtaacacgccacatcttgcgaatatatgtgtagaaactgccggaaatcgtcgtggtattcactccagagcgatgaaaacgtttcagtttgctcatggaaaacggtgtaacaagggtgaacactatcccatatcaccagctcaccgtctttcattgccatacggaactccggatgagcattcatcaggcgggcaagaatgtgaataaaggccggataaaacttgtgcttatttttctttacggtctttaaaaaggccgtaatatccagctgaacggtctggttataggtacattgagcaactgactgaaatgcctcaaaatgttctttacgatgccattgggatatatcaacggtggtatatccagtgatttttttctccat")),

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
