//! A library of known sequences we can use to automatically add features to a sequence.
//! This will eventually load data from an online source. Currently, it includes some common items,
//! and is stored in memory.
//!
//! Some common genes: https://en.vectorbuilder.com/search/gene.html

use na_seq::{
    Nucleotide, Seq,
    amino_acids::{AminoAcid, CodingResult},
    seq_complement, seq_from_str,
};

use crate::{
    misc_types::{
        Feature, FeatureDirection,
        FeatureType::{
            self, AntibioticResistance, CodingRegion, Generic, Ori, Promoter, ProteinBind,
            RibosomeBindSite, Terminator,
        },
    },
    reading_frame::ReadingFrame,
    util::{RangeIncl, match_subseq},
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
        let offset = orf.offset();
        let seq_ = orf.arrange_seq(seq);
        let len = seq_.len();

        let mut tag_open = None; // Inner: Start index.

        for i_ in 0..len / 3 {
            let i = i_ * 3; // The actual sequence index.

            let nts = &seq_[i..i + 3];

            let mut matched = false;
            if let CodingResult::AminoAcid(aa) = AminoAcid::from_codons(nts.try_into().unwrap()) {
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
                            // RangeIncl::new(tag_open.unwrap() + 1 + offset, i + 2 + offset)
                            RangeIncl::new(tag_open.unwrap() + 1 + offset, i + offset)
                        }
                        _ => RangeIncl::new(
                            seq_len_full - (i + 2 + offset),
                            seq_len_full - (tag_open.unwrap() + offset) - 1,
                        ),
                    };

                    // todo: Likely you want the appropriate direction, determined by the reading frame.
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

// todo: Map Aa sequences in addition to DNA seqs; more general.

/// Find common promoters and Oris.
/// The order of the written sequences matters: It determines direction.
fn find_misc(seq: &[Nucleotide]) -> Vec<Feature> {
    let items = vec![
        FeatureMapItem::new(
            "AmpR promoter",
            Promoter,
            seq_from_str(
                "cgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataacattgaaaaaggaagagt",
            ),
        ),
        // Perhaps slightly different from the above AmpR promoter?
        FeatureMapItem::new(
            "AmpR promoter",
            Promoter,
            seq_from_str(
                "cgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagt",
            ),
        ),
        FeatureMapItem::new("T3 promoter", Promoter, seq_from_str("aattaaccctcactaaagg")),
        FeatureMapItem::new(
            "lac promoter",
            Promoter,
            seq_from_str("tttacactttatgcttccggctcgtatgttg"),
        ),
        FeatureMapItem::new(
            "AmpR promoter",
            Promoter,
            seq_from_str(
                "actcttcctttttcaatattattgaagcatttatcagggttattgtctcatgagcggatacatatttgaatgtatttagaaaaataaacaaataggggttccgcg",
            ),
        ),
        FeatureMapItem::new(
            "ARA promoter",
            Promoter,
            seq_from_str("ctgacgctttttatcgcaactctctactg"),
        ),
        FeatureMapItem::new(
            "lac promoter",
            Promoter,
            seq_from_str("CAACATACGAGCCGGAAGCATAAAGTGTAAA"),
        ),
        FeatureMapItem::new(
            "cat promoter",
            Promoter,
            seq_from_str(
                "TGATCGGCACGTAAGAGGKTCCAACTTTCACCATAATGAAATAAGATCACTACCGGGCGTATTTTTTGAGTTRTCGAGATTTTCAGGAGCTAAGGAAGCTAAA",
            ),
        ),
        FeatureMapItem::new(
            "lacI promoter",
            Promoter,
            seq_from_str(
                "GACACCATCGAATGGCGCAAAACCTTTCGCGGTATGGCATGATAGCGCCCGGAAGAGAGTCAATTCAGGGTGGTGAAT",
            ),
        ),
        FeatureMapItem::new(
            "lacIq promoter",
            Promoter,
            seq_from_str(
                "gacaccatcgaatggtgcaaaacctttcgcggtatggcatgatagcgcccggaagagagtcaattcagggtggtgaat",
            ),
        ),
        FeatureMapItem::new(
            "tac promoter",
            Promoter,
            seq_from_str("TTGACAATTAATCATCGGCTCGTATAATG"),
        ),
        FeatureMapItem::new(
            "trc promoter",
            Promoter,
            seq_from_str("TTGACAATTAATCATCCGGCYCGTATAATG"),
        ),
        FeatureMapItem::new(
            "CMV promoter",
            Promoter,
            seq_from_str(
                "gtgatgcggttttggcagtacatcaatgggcgtggatagcggtttgactcacggggatttccaagtctccaccccattgacgtcaatgggagtttgttttggcaccaaaatcaacgggactttccaaaatgtcgtaacaactccgccccattgacgcaaatgggcggtaggcgtgtacggtgggaggtctatataagcagagct",
            ),
        ),
        FeatureMapItem::new(
            "T37 promoter",
            Promoter,
            seq_from_str("aattaaccctcactaaagg"),
        ),
        FeatureMapItem::new("T7 promoter", Promoter, seq_from_str("taatacgactcactatagg")),
        FeatureMapItem::new(
            "Cat promoter",
            Promoter,
            seq_complement(&seq_from_str(
                "tttagcttccttagctcctgaaaatctcgataactcaaaaaatacgcccggtagtgatcttatttcattatggtgaaagttggaacctcttacgtgccgatca",
            )),
        ),
        FeatureMapItem::new(
            "sacB promoter",
            Promoter,
            seq_from_str(
                "cacatatacctgccgttcactattatttagtgaaatgagatattatgatattttctgaattgtgattaaaaaggcaactttatgcccatgcaacagaaactataaaaaatacagagaatgaaaagaaacagatagattttttagttctttaggcccgtagtctgcaaatccttttatgattttctatcaaacaaaagaggaaaatagaccagttgcaatccaaacgagagtctaatagaatgaggtcgaaaagtaaatcgcgcgggtttgttactgataaagcaggcaagacctaaaatgtgtaaagggcaaagtgtatactttggcgtcaccccttacatattttaggtctttttttattgtgcgtaactaacttgccatcttcaaacaggagggctggaagaagcagaccgctaacacagtacataaaaaaggagacatgaacg",
            ),
        ),
        FeatureMapItem::new(
            "lacZα",
            CodingRegion,
            seq_complement(&seq_from_str(
                "CTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGTAATCATGGTCAT",
            )),
        ),
        FeatureMapItem::new(
            "lacZα",
            CodingRegion,
            seq_from_str(
                "ctatgcggcatcagagcagattgtactgagagtgcaccatatgcggtgtgaaataccgcacagatgcgtaaggagaaaataccgcatcaggcgccattcgccattcaggctgcgcaactgttgggaagggcgatcggtgcgggcctcttcgctattacgccagctggcgaaagggggatgtgctgcaaggcgattaagttgggtaacgccagggttttcccagtcacgacgttgtaaaacgacggccagtgaattcgagctcggtacccggggatcctctagagtcgacctgcaggcatgcaagcttggcgtaatcatggtcat",
            ),
        ),
        FeatureMapItem::new(
            "lacZα",
            CodingRegion,
            seq_complement(&seq_from_str(
                "ctaatcaagttttttggggtcgaggtgccgtaaagcactaaatcggaaccctaaagggagcccccgatttagagcttgacggggaaagccggcgaacgtggcgagaaaggaagggaagaaagcgaaaggagcgggcgctagggcgctggcaagtgtagcggtcacgctgcgcgtaaccaccacacccgccgcgcttaatgcgccgctacagggcgcgtcccattcgccattcaggctgcgcaactgttgggaagggcgatcggtgcgggcctcttcgctattacgccagctggcgaaagggggatgtgctgcaaggcgattaagttgggtaacgccagggttttcccagtcacgacgttgtaaaacgacggccagtgagcgcgcgtaatacgactcactatagggcgaattgggtaccgggccccccctcgaggtcgacggtatcgataagcttgatatcgaattcctgcagcccgggggatccactagttctagagcggccgccaccgcggtggagctccagcttttgttccctttagtgagggttaattgcgcgcttggcgtaatcatggtcat",
            )),
        ),
        FeatureMapItem::new(
            "lacZα",
            CodingRegion,
            seq_from_str(
                "atgaccatgattacggattcactggccgtcgttttacaacgtcgtgactgggaaaaccctggcgttacccaacttaatcgccttgcagcacatccccctttcgccagctggcgtaatagcgaagaggcccgcaccgatcgcccttcccaacagttgcgcagcctgaatggcgaatggcgctttgcctggtttccggcaccagaagcggtgccggaaagctggctggagtgcgatcttcctgaggccgatactgtcgtcgtcccctcaaactggcagatgcacggttacgatgcgcccatctacaccaacgtaacctatcccattacggtcaatccgccgtttgttcccacggagaatccgacgggttgttactcgctcacatttaatgttgatgaaagctggctacaggaaggccagacgcgaattatttttgatggcgttggaatt",
            ),
        ),
        FeatureMapItem::new(
            "lacI",
            CodingRegion,
            seq_complement(&seq_from_str(
                "tcactgcccgctttccagtcgggaaacctgtcgtgccagctgcattaatgaatcggccaacgcgcggggagaggcggtttgcgtattgggcgccagggtggtttttcttttcaccagtgagacgggcaacagctgattgcccttcaccgcctggccctgagagagttgcagcaagcggtccacgctggtttgccccagcaggcgaaaatcctgtttgatggtggttaacggcgggatataacatgagctgtcttcggtatcgtcgtatcccactaccgagatgtccgcaccaacgcgcagcccggactcggtaatggcgcgcattgcgcccagcgccatctgatcgttggcaaccagcatcgcagtgggaacgatgccctcattcagcatttgcatggtttgttgaaaaccggacatggcactccagtcgccttcccgttccgctatcggctgaatttgattgcgagtgagatatttatgccagccagccagacgcagacgcgccgagacagaacttaatgggcccgctaacagcgcgatttgctggtgacccaatgcgaccagatgctccacgcccagtcgcgtaccgtcttcatgggagaaaataatactgttgatgggtgtctggtcagagacatcaagaaataacgccggaacattagtgcaggcagcttccacagcaatggcatcctggtcatccagcggatagttaatgatcagcccactgacgcgttgcgcgagaagattgtgcaccgccgctttacaggcttcgacgccgcttcgttctaccatcgacaccaccacgctggcacccagttgatcggcgcgagatttaatcgccgcgacaatttgcgacggcgcgtgcagggccagactggaggtggcaacgccaatcagcaacgactgtttgcccgccagttgttgtgccacgcggttgggaatgtaattcagctccgccatcgccgcttccactttttcccgcgttttcgcagaaacgtggctggcctggttcaccacgcgggaaacggtctgataagagacaccggcatactctgcgacatcgtataacgttactggtttcac",
            )),
        ),
        FeatureMapItem::new(
            "lacI",
            CodingRegion,
            seq_from_str(
                "gtgaaaccagtaacgttatacgatgtcgcagagtatgccggtgtctcttatcagaccgtttcccgcgtggtgaaccaggccagccacgtttctgcgaaaacgcgggaaaaagtggaagcggcgatggcggagctgaattacattcccaaccgcgtggcacaacaactggcgggcaaacagtcgttgctgattggcgttgccacctccagtctggccctgcacgcgccgtcgcaaattgtcgcggcgattaaatctcgcgccgatcaactgggtgccagcgtggtggtgtcgatggtagaacgaagcggcgtcgaagcctgtaaagcggcggtgcacaatcttctcgcgcaacgcgtcagtgggctgatcattaactatccgctggatgaccaggatgccattgctgtggaagctgcctgcactaatgttccggcgttatttcttgatgtctctgaccagacacccatcaacagtattattttctcccatgaagacggtacgcgactgggcgtggagcatctggtcgcattgggtcaccagcaaatcgcgctgttagcgggcccattaagttctgtctcggcgcgtctgcgtctggctggctggcataaatatctcactcgcaatcaaattcagccgatagcggaacgggaaggcgactggagtgccatgtccggttttcaacaaaccatgcaaatgctgaatgagggcatcgttcccactgcgatgctggttgccaacgatcagatggcgctgggcgcaatgcgcgccattaccgagtccgggctgcgcgttggtgcggatatctcggtagtgggatacgacgataccgaagacagctcatgttatatcccgccgttaaccaccatcaaacaggattttcgcctgctggggcaaaccagcgtggaccgcttgctgcaactctctcagggccaggcggtgaagggcaatcagctgttgcccgtctcactggtgaaaagaaaaaccaccctggcgcccaatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctggcacgacaggtttcccgactggaaagcgggcagtga",
            ),
        ),
        FeatureMapItem::new(
            "lacI",
            CodingRegion,
            seq_from_str(
                "gtgaaaccagtaacgttatacgatgtcgcagagtatgccggtgtctcttatcagaccgtttcccgcgtggtgaaccaggccagccacgtttctgcgaaaacgcgggaaaaagtggaagcggcgatggcggagctgaattacattcccaaccgcgtggcacaacaactggcgggcaaacagtcgttgctgattggcgttgccacctccagtctggccctgcacgcgccgtcgcaaattgtcgcggcgattaaatctcgcgccgatcaactgggtgccagcgtggtggtgtcgatggtagaacgaagcggcgtcgaagcctgtaaagcggcggtgcacaatcttctcgcgcaacgcgtcagtgggctgatcattaactatccgctggatgaccaggatgccattgctgtggaagctgcctgcactaatgttccggcgttatttcttgatgtctctgaccagacacccatcaacagtattattttctcccatgaagacggtacgcgactgggcgtggagcatctggtcgcattgggtcaccagcaaatcgcgctgttagcgggcccattaagttctgtctcggcgcgtctgcgtctggctggctggcataaatatctcactcgcaatcaaattcagccgatagcggaacgggaaggcgactggagtgccatgtccggttttcaacaaaccatgcaaatgctgaatgagggcatcgttcccactgcgatgctggttgccaacgatcagatggcgctgggcgcaatgcgcgccattaccgagtccgggctgcgcgttggtgcggatatctcggtagtgggatacgacgataccgaagacagctcatgttatatcccgccgtcaaccaccatcaaacaggattttcgcctgctggggcaaaccagcgtggaccgcttgctgcaactctctcagggccaggcggtgaagggcaatcagctgttgcccgtctcactggtgaaaagaaaaaccaccctggcgcccaatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctggcacgacaggtttcccgactggaaagcgggcagtga",
            ),
        ),
        FeatureMapItem::new(
            "araC",
            CodingRegion,
            seq_complement(&seq_from_str(
                "ttatgacaacttgacggctacatcattcactttttcttcacaaccggcacggaactcgctcgggctggccccggtgcattttttaaatacccgcgagaaatagagttgatcgtcaaaaccaacattgcgaccgacggtggcgataggcatccgggtggtgctcaaaagcagcttcgcctggctgatacgttggtcctcgcgccagcttaagacgctaatccctaactgctggcggaaaagatgtgacagacgcgacggcgacaagcaaacatgctgtgcgacgctggcgatatcaaaattgctgtctgccaggtgatcgctgatgtactgacaagcctcgcgtacccgattatccatcggtggatggagcgactcgttaatcgcttccatgcgccgcagtaacaattgctcaagcagatttatcgccagcagctccgaatagcgcccttccccttgcccggcgttaatgatttgcccaaacaggtcgctgaaatgcggctggtgcgcttcatccgggcgaaagaaccccgtattggcaaatattgacggccagttaagccattcatgccagtaggcgcgcggacgaaagtaaacccactggtgataccattcgcgagcctccggatgacgaccgtagtgatgaatctctcctggcgggaacagcaaaatatcacccggtcggcaaacaaattctcgtccctgatttttcaccaccccctgaccgcgaatggtgagattgagaatataacctttcattcccagcggtcggtcgataaaaaaatcgagataaccgttggcctcaatcggcgttaaacccgccaccagatgggcattaaacgagtatcccggcagcaggggatcattttgcgcttcagccat",
            )),
        ),
        FeatureMapItem::new(
            "lac operator",
            ProteinBind,
            seq_from_str("ttgtgagcggataacaa"),
        ),
        FeatureMapItem::new(
            "lac operator",
            ProteinBind,
            seq_from_str("ttgttatccgctcacaa"),
        ),
        FeatureMapItem::new(
            "RBS",
            RibosomeBindSite,
            seq_from_str("tttgtttaactttaagaaggaga"),
        ),
        FeatureMapItem::new(
            "T7 terminator",
            Terminator,
            seq_from_str("ctagcataaccccttggggcctctaaacgggtcttgaggggttttttg"),
        ),
        FeatureMapItem::new(
            "rrnB T1 terminator",
            Terminator,
            seq_from_str(
                "caaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctcctgagtaggacaaat",
            ),
        ),
        FeatureMapItem::new(
            "rrnB T2 terminator",
            Terminator,
            seq_from_str("agaaggccatcctgacggatggcctttt"),
        ),
        FeatureMapItem::new(
            "SV40 promoter",
            Promoter,
            seq_from_str(
                "ctgaggcggaaagaaccagctgtggaatgtgtgtcagttagggtgtggaaagtccccaggctccccagcaggcagaagtatgcaaagcatgcatctcaattagtcagcaaccaggtgtggaaagtccccaggctccccagcaggcagaagtatgcaaagcatgcatctcaattagtcagcaaccatagtcccgcccctaactccgcccatcccgcccctaactccgcccagttccgcccattctccgccccatggctgactaattttttttatttatgcagaggccgaggccgcctcggcctctgagctattccagaagtagtgaggaggcttttttggaggcctaggcttttgcaaa",
            ),
        ),
        FeatureMapItem::new(
            "CAP binding site",
            ProteinBind,
            seq_from_str("atgagtgagctaactcacatta"),
        ),
        FeatureMapItem::new(
            "CAP binding site",
            ProteinBind,
            seq_from_str("atgagtgagctaactcacatta"),
        ),
        FeatureMapItem::new(
            "pHluorin2",
            CodingRegion,
            seq_from_str(
                "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGAGCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACGAGCACCTGGTGTACATCATGGCCGACAAGCAGAAGAACGGCACCAAGGCCATCTTCCAGGTGCACCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGCACACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCACGGCATGGACGAGCTGTACAAGTAA",
            ),
        ),
        FeatureMapItem::new(
            "SacB",
            CodingRegion,
            seq_from_str(
                "atgaacatcaaaaagtttgcaaaacaagcaacagtattaacctttactaccgcactgctggcaggaggcgcaactcaagcgtttgcgaaagaaacgaaccaaaagccatataaggaaacatacggcatttcccatattacacgccatgatatgctgcaaatccctgaacagcaaaaaaatgaaaaatataaagttcctgagttcgattcgtccacaattaaaaatatctcttctgcaaaaggcctggacgtttgggacagctggccattacaaaacactgacggcactgtcgcaaactatcacggctaccacatcgtctttgcattagccggagatcctaaaaatgcggatgacacatcgatttacatgttctatcaaaaagtcggcgaaacttctattgacagctggaaaaacgctggccgcgtctttaaagacagcgacaaattcgatgcaaatgattctatcctaaaagaccaaacacaagaatggtcaggttcagccacatttacatctgacggaaaaatccgtttattctacactgatttctccggtaaacattacggcaaacaaacactgacaactgcacaagttaacgtatcagcatcagacagctctttgaacatcaacggtgtagaggattataaatcaatctttgacggtgacggaaaaacgtatcaaaatgtacagcagttcatcgatgaaggcaactacagctcaggcgacaaccatacgctgagagatcctcactacgtagaagataaaggccacaaatacttagtatttgaagcaaacactggaactgaagatggctaccaaggcgaagaatctttatttaacaaagcatactatggcaaaagcacatcattcttccgtcaagaaagtcaaaaacttctgcaaagcgataaaaaacgcacggctgagttagcaaacggcgctctcggtatgattgagctaaacgatgattacacactgaaaaaagtgatgaaaccgctgattgcatctaacacagtaacagatgaaattgaacgcgcgaacgtctttaaaatgaacggcaaatggtacctgttcactgactcccgcggatcaaaaatgacgattgacggcattacgtctaacgatatttacatgcttggttatgtttctaattctttaactggcccatacaagccgctgaacaaaactggccttgtgttaaaaatggatcttgatcctaacgatgtaacctttacttactcacacttcgctgtacctcaagcgaaaggaaacaatgtcgtgattacaagctatatgacaaacagaggattctacgcagacaaacaatcaacgtttgcgcctagcttcctgctgaacatcaaaggcaagaaaacatctgttgtcaaagacagcatccttgaacaaggacaattaacagttaacaaataa",
            ),
        ),
        FeatureMapItem::new(
            "AmpR",
            AntibioticResistance,
            seq_complement(&seq_from_str(
                "gttaccaatgcttaatcagtgaggcacctatctcagcgatctgtctatttcgttcatccatagttgcctgactccccgtcgtgtagataactacgatacgggagggcttaccatctggccccagtgctgcaatgataccgcgagacccacgctcaccggctccagatttatcagcaataaaccagccagccggaagggccgagcgcagaagtggtcctgcaactttatccgcctccatccagtctattaattgttgccgggaagctagagtaagtagttcgccagttaatagtttgcgcaacgttgttgccattgctacaggcatcgtggtgtcacgctcgtcgtttggtatggcttcattcagctccggttcccaacgatcaaggcgagttacatgatcccccatgttgtgcaaaaaagcggttagctccttcggtcctccgatcgttgtcagaagtaagttggccgcagtgttatcactcatggttatggcagcactgcataattctcttactgtcatgccatccgtaagatgcttttctgtgactggtgagtactcaaccaagtcattctgagaatagtgtatgcggcgaccgagttgctcttgcccggcgtcaatacgggataataccgcgccacatagcagaactttaaaagtgctcatcattggaaaacgttcttcggggcgaaaactctcaaggatcttaccgctgttgagatccagttcgatgtaacccactcgtgcacccaactgatcttcagcatcttttactttcaccagcgtttctgggtgagcaaaaacaggaaggcaaaatgccgcaaaaaagggaataagggcgacacggaaatgttgaatactcat",
            )),
        ),
        FeatureMapItem::new(
            "AmpR",
            AntibioticResistance,
            seq_from_str(
                "atgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgtagcaatggcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtcgcgcggtatcattgcagcactggggccagatggtaagccctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaa",
            ),
        ),
        FeatureMapItem::new(
            "AmpR",
            AntibioticResistance,
            seq_from_str(
                "atgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtgttgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgcagcaatggcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagccctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaa",
            ),
        ),
        FeatureMapItem::new(
            "AmpR",
            AntibioticResistance,
            seq_complement(&seq_from_str(
                "ttaccaatgcttaatcagtgaggcacctatctcagcgatctgtctatttcgttcatccatagttgcctgactccccgtcgtgtagataactacgatacgggagggcttaccatctggccccagtgctgcaatgataccgcgagacccacgctcaccggctccagatttatcagcaataaaccagccagccggaagggccgagcgcagaagtggtcctgcaactttatccgcctccatccagtctattaattgttgccgggaagctagagtaagtagttcgccagttaatagtttgcgcaacgttgttgccattgctgcaggcatcgtggtgtcacgctcgtcgtttggtatggcttcattcagctccggttcccaacgatcaaggcgagttacatgatcccccatgttgtgcaaaaaagcggttagctccttcggtcctccgatcgttgtcagaagtaagttggccgcagtgttatcactcatggttatggcagcactgcataattctcttactgtcatgccatccgtaagatgcttttctgtgactggtgagtactcaaccaagtcattctgagaatagtgtatgcggcgaccgagttgctcttgcccggcgtcaatacgggataataccgcgccacatagcagaactttaaaagtgctcatcattggaaaacgttcttcggggcgaaaactctcaaggatcttaccgctgttgagatccagttcgatgtaacccactcgtgcacccaactgatcttcagcatcttttactttcaccagcgtttctgggtg",
            )),
        ),
        FeatureMapItem::new(
            "KanR",
            AntibioticResistance,
            seq_from_str(
                "gttagaaaaactcatcgagcatcaaatgaaactgcaatttattcatatcaggattatcaataccatatttttgaaaaagccgtttctgtaatgaaggagaaaactcaccgaggcagttccataggatggcaagatcctggtatcggtctgcgattccgactcgtccaacatcaatacaacctattaatttcccctcgtcaaaaataaggttatcaagtgagaaatcaccatgagtgacgactgaatccggtgagaatggcaaaagtttatgcatttctttccagacttgttcaacaggccagccattacgctcgtcatcaaaatcactcgcatcaaccaaaccgttattcattcgtgattgcgcctgagcgagacgaaatacgcgatcgctgttaaaaggacaattacaaacaggaatcgaatgcaaccggcgcaggaacactgccagcgcatcaacaatattttcacctgaatcaggatattcttctaatacctggaatgctgttttcccagggatcgcagtggtgagtaaccatgcatcatcaggagtacggataaaatgcttgatggtcggaagaggcataaattccgtcagccagtttagtctgaccatctcatctgtaacatcattggcaacgctacctttgccatgtttcagaaacaactctggcgcatcgggcttcccatacaatcgatagattgtcgcacctgattgcccgacattatcgcgagcccatttatacccatataaatcagcatccatgttggaatttaatcgcggcctagagcaagacgtttcccgttgaatatggctcat",
            ),
        ),
        FeatureMapItem::new(
            "KanR",
            AntibioticResistance,
            seq_from_str(
                "atgagccatattcaacgggaaacgtcttgctctaggccgcgattaaattccaacatggatgctgatttatatgggtataaatgggctcgcgataatgtcgggcaatcaggtgcgacaatctatcgattgtatgggaagcccgatgcgccagagttgtttctgaaacatggcaaaggtagcgttgccaatgatgttacagatgagatggtcagactaaactggctgacggaatttatgcctcttccgaccatcaagcattttatccgtactcctgatgatgcatggttactcaccactgcgatccccgggaaaacagcattccaggtattagaagaatatcctgattcaggtgaaaatattgttgatgcgctggcagtgttcctgcgccggttgcattcgattcctgtttgtaattgtccttttaacagcgatcgcgtatttcgtctcgctcaggcgcaatcacgaatgaataacggtttggttgatgcgagtgattttgatgacgagcgtaatggctggcctgttgaacaagtctggaaagaaatgcataaacttttgccattctcaccggattcagtcgtcactcatggtgatttctcacttgataaccttatttttgacgaggggaaattaataggttgtattgatgttggacgagtcggaatcgcagaccgataccaggatcttgccatcctatggaactgcctcggtgagttttctccttcattacagaaacggctttttcaaaaatatggtattgataatcctgatatgaataaattgcagtttcatttgatgctcgatgagtttttctaa",
            ),
        ),
        FeatureMapItem::new(
            "CmR",
            AntibioticResistance,
            seq_complement(&seq_from_str(
                "attacgccccgccctgccactcatcgcagtactgttgtaattcattaagcattctgccgacatggaagccatcacaaacggcatgatgaacctgaatcgccagcggcatcagcaccttgtcgccttgcgtataatatttgcccatagtgaaaacgggggcgaagaagttgtccatattggccacgtttaaatcaaaactggtgaaactcacccagggattggctgagacgaaaaacatattctcaataaaccctttagggaaataggccaggttttcaccgtaacacgccacatcttgcgaatatatgtgtagaaactgccggaaatcgtcgtggtattcactccagagcgatgaaaacgtttcagtttgctcatggaaaacggtgtaacaagggtgaacactatcccatatcaccagctcaccgtctttcattgccatacggaactccggatgagcattcatcaggcgggcaagaatgtgaataaaggccggataaaacttgtgcttatttttctttacggtctttaaaaaggccgtaatatccagctgaacggtctggttataggtacattgagcaactgactgaaatgcctcaaaatgttctttacgatgccattgggatatatcaacggtggtatatccagtgatttttttctccat",
            )),
        ),
        FeatureMapItem::new(
            "CmR",
            AntibioticResistance,
            seq_complement(&seq_from_str(
                "ttacgccccgccctgccactcatcgcagtactgttgtaattcattaagcattctgccgacatggaagccatcacaaacggcatgatgaacctgaatcgccagcggcatcagcaccttgtcgccttgcgtataatatttgcccatagtgaaaacgggggcgaagaagttgtccatattggccacgtttaaatcaaaactggtgaaactcacccagggattggctgagacgaaaaacatattctcaataaaccctttagggaaataggccaggttttcaccgtaacacgccacatcttgcgaatatatgtgtagaaactgccggaaatcgtcgtggtattcactccagagcgatgaaaacgtttcagtttgctcatggaaaacggtgtaacaagggtgaacactatcccatatcaccagctcaccgtctttcattgccatacggaactccggatgagcattcatcaggcgggcaagaatgtgaataaaggccggataaaacttgtgcttatttttctttacggtctttaaaaaggccgtaatatccagctgaacggtctggttataggtacattgagcaactgactgaaatgcctcaaaatgttctttacgatgccattgggatatatcaacggtggtatatccagtgatttttttctccat",
            )),
        ),
        FeatureMapItem::new(
            "NeoR/KanR",
            AntibioticResistance,
            seq_from_str(
                "atgattgaacaagatggattgcacgcaggttctccggccgcttgggtggagaggctattcggctatgactgggcacaacagacaatcggctgctctgatgccgccgtgttccggctgtcagcgcaggggcgcccggttctttttgtcaagaccgacctgtccggtgccctgaatgaactgcaagacgaggcagcgcggctatcgtggctggccacgacgggcgttccttgcgcagctgtgctcgacgttgtcactgaagcgggaagggactggctgctattgggcgaagtgccggggcaggatctcctgtcatctcaccttgctcctgccgagaaagtatccatcatggctgatgcaatgcggcggctgcatacgcttgatccggctacctgcccattcgaccaccaagcgaaacatcgcatcgagcgagcacgtactcggatggaagccggtcttgtcgatcaggatgatctggacgaagagcatcaggggctcgcgccagccgaactgttcgccaggctcaaggcgagcatgcccgacggcgaggatctcgtcgtgacccatggcgatgcctgcttgccgaatatcatggtggaaaatggccgcttttctggattcatcgactgtggccggctgggtgtggcggaccgctatcaggacatagcgttggctacccgtgatattgctgaagagcttggcggcgaatgggctgaccgcttcctcgtgctttacggtatcgccgctcccgattcgcagcgcatcgccttctatcgccttcttgacgagttcttctga",
            ),
        ),
        FeatureMapItem::new(
            "CMV enhancer",
            CodingRegion,
            seq_from_str(
                "cgttacataacttacggtaaatggcccgcctggctgaccgcccaacgacccccgcccattgacgtcaataatgacgtatgttcccatagtaacgccaatagggactttccattgacgtcaatgggtggagtatttacggtaaactgcccacttggcagtacatcaagtgtatcatatgccaagtacgccccctattgacgtcaatgacggtaaatggcccgcctggcattatgcccagtacatgaccttatgggactttcctacttggcagtacatctacgtattagtcatcgctattaccatg",
            ),
        ),
        FeatureMapItem::new(
            "EGFP",
            CodingRegion,
            seq_from_str(
                "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA",
            ),
        ),
        FeatureMapItem::new(
            "TcR",
            CodingRegion,
            seq_from_str(
                "atgaaatctaacaatgcgctcatcgtcatcctcggcaccgtcaccctggatgctgtaggcataggcttggttatgccggtactgccgggcctcttgcgggatatcgtccattccgacagcatcgccagtcactatggcgtgctgctagcgctatatgcgttgatgcaatttctatgcgcacccgttctcggagcactgtccgaccgctttggccgccgcccagtcctgctcgcttcgctacttggagccactatcgactacgcgatcatggcgaccacacccgtcctgtggatcctctacgccggacgcatcgtggccggcatcaccggcgccacaggtgcggttgctggcgcctatatcgccgacatcaccgatggggaagatcgggctcgccacttcgggctcatgagcgcttgtttcggcgtgggtatggtggcaggccccgtggccgggggactgttgggcgccatctccttgcatgcaccattccttgcggcggcggtgctcaacggcctcaacctactactgggctgcttcctaatgcaggagtcgcataagggagagcgtcgaccgatgcccttgagagccttcaacccagtcagctccttccggtgggcgcggggcatgactatcgtcgccgcacttatgactgtcttctttatcatgcaactcgtaggacaggtgccggcagcgctctgggtcattttcggcgaggaccgctttcgctggagcgcgacgatgatcggcctgtcgcttgcggtattcggaatcttgcacgccctcgctcaagccttcgtcactggtcccgccaccaaacgtttcggcgagaagcaggccattatcgccggcatggcggccgacgcgctgggctacgtcttgctggcgttcgcgacgcgaggctggatggccttccccattatgattcttctcgcttccggcggcatcgggatgcccgcgttgcaggccatgctgtccaggcaggtagatgacgaccatcagggacagcttcaaggatcgctcgcggctcttaccagcctaacttcgatcattggaccgctgatcgtcacggcgatttatgccgcctcggcgagcacatggaacgggttggcatggattgtaggcgccgccctataccttgtctgcctccccgcgttgcgtcgcggtgcatggagccgggccacctcgacctga",
            ),
        ),
        FeatureMapItem::new(
            "GST",
            CodingRegion,
            seq_from_str(
                "atgtcccctatactaggttattggaaaattaagggccttgtgcaacccactcgacttcttttggaatatcttgaagaaaaatatgaagagcatttgtatgagcgcgatgaaggtgataaatggcgaaacaaaaagtttgaattgggtttggagtttcccaatcttccttattatattgatggtgatgttaaattaacacagtctatggccatcatacgttatatagctgacaagcacaacatgttgggtggttgtccaaaagagcgtgcagagatttcaatgcttgaaggagcggttttggatattagatacggtgtttcgagaattgcatatagtaaagactttgaaactctcaaagttgattttcttagcaagctacctgaaatgctgaaaatgttcgaagatcgtttatgtcataaaacatatttaaatggtgatcatgtaacccatcctgacttcatgttgtatgacgctcttgatgttgttttatacatggacccaatgtgcctggatgcgttcccaaaattagtttgttttaaaaaacgtattgaagctatcccacaaattgataagtacttgaaatccagcaagtatatagcatggcctttgcagggctggcaagccacgtttggtggtggcgaccatcctccaaaa",
            ),
        ),
        // todo: How general is this?
        FeatureMapItem::new(
            "Ori",
            Ori,
            seq_complement(&seq_from_str(
                "ttttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaagaacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaa",
            )),
        ),
        FeatureMapItem::new(
            "Ori",
            Ori,
            seq_complement(&seq_from_str(
                "tttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaaggacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaa",
            )),
        ),
        FeatureMapItem::new(
            "f1 ori",
            Ori,
            seq_from_str(
                "acgcgccctgtagcggcgcattaagcgcggcgggtgtggtggttacgcgcagcgtgaccgctacacttgccagcgccctagcgcccgctcctttcgctttcttcccttcctttctcgccacgttcgccggctttccccgtcaagctctaaatcgggggctccctttagggttccgatttagtgctttacggcacctcgaccccaaaaaacttgattagggtgatggttcacgtagtgggccatcgccctgatagacggtttttcgccctttgacgttggagtccacgttctttaatagtggactcttgttccaaactggaacaacactcaaccctatctcggtctattcttttgatttataagggattttgccgatttcggcctattggttaaaaaatgagctgatttaacaaaaatttaacgcgaattttaacaaaatattaacgcttacaattt",
            ),
        ),
        FeatureMapItem::new(
            "f1 ori",
            Ori,
            seq_complement(&seq_from_str(
                "aaattgtaaacgttaatattttgttaaaattcgcgttaaatttttgttaaatcagctcattttttaaccaataggccgaaatcggcaaaatcccttataaatcaaaagaatagaccgagatagggttgagtgttgttccagtttggaacaagagtccactattaaagaacgtggactccaacgtcaaagggcgaaaaaccgtctatcagggcgatggcccactacgtgaaccatcaccctaatcaagttttttggggtcgaggtgccgtaaagcactaaatcggaaccctaaagggagcccccgatttagagcttgacggggaaagccggcgaacgtggcgagaaaggaagggaagaaagcgaaaggagcgggcgctagggcgctggcaagtgtagcggtcacgctgcgcgtaaccaccacacccgccgcgcttaatgcgccgctacagggcgcgt",
            )),
        ),
        FeatureMapItem::new(
            "f1 ori",
            Ori,
            seq_from_str(
                "acgcgccctgtagcggcgcattaagcgcggcgggtgtggtggttacgcgcagcgtgaccgctacacttgccagcgccctagcgcccgctcctttcgctttcttcccttcctttctcgccacgttcgccggctttccccgtcaagctctaaatcgggggctccctttagggttccgatttagtgctttacggcacctcgaccccaaaaaacttgatttgggtgatggttcacgtagtgggccatcgccctgatagacggtttttcgccctttgacgttggagtccacgttctttaatagtggactcttgttccaaactggaacaacactcaaccctatctcgggctattcttttgatttataagggattttgccgatttcggcctattggttaaaaaatgagctgatttaacaaaaatttaacgcgaattttaacaaaatattaacgtttacaattt",
            ),
        ),
        FeatureMapItem::new(
            "P15A ori",
            Ori,
            seq_complement(&seq_from_str(
                "ttttccataggctccgcccccctgacaagcatcacgaaatctgacgctcaaatcagtggtggcgaaacccgacaggactataaagataccaggcgtttccccctggcggctccctcgtgcgctctcctgttcctgcctttcggtttaccggtgtcattccgctgttatggccgcgtttgtctcattccacgcctgacactcagttccgggtaggcagttcgctccaagctggactgtatgcacgaaccccccgttcagtccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggaaagacatgcaaaagcaccactggcagcagccactggtaattgatttagaggagttagtcttgaagtcatgcgccggttaaggctaaactgaaaggacaagttttggtgactgcgctcctccaagccagttacctcggttcaaagagttggtagctcagagaaccttcgaaaaaccgccctgcaaggcggttttttcgttttcagagcaagagattacgcgcagaccaaaacgatctcaa",
            )),
        ),
        FeatureMapItem::new(
            "SV40 ori",
            Ori,
            seq_from_str(
                "atcccgcccctaactccgcccagttccgcccattctccgccccatggctgactaattttttttatttatgcagaggccgaggccgcctcggcctctgagctattccagaagtagtgaggaggcttttttggaggcc",
            ),
        ),
        FeatureMapItem::new(
            "ROP",
            CodingRegion,
            seq_complement(&seq_from_str(
                "tcagaggttttcaccgtcatcaccgaaacgcgcgaggcagctgcggtaaagctcatcagcgtggtcgtgaagcgattcacagatgtctgcctgttcatccgcgtccagctcgttgagtttctccagaagcgttaatgtctggcttctgataaagcgggccatgttaagggcggttttttcctgtttggtcac",
            )),
        ),
        FeatureMapItem::new(
            "BOM",
            Generic,
            seq_from_str(
                "cctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcaatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagtatacactccgctatcgctacgtgactgggtcatggctgcg",
            ),
        ),
        FeatureMapItem::new(
            "BOM",
            Generic,
            seq_from_str(
                "cctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcactggtgcactctcagtacaatctgctctgatgccgcatagttaagccagtatacactccgctatcgctacgtgactgggtcatggctgcg",
            ),
        ),
        FeatureMapItem::new(
            "BOM",
            Generic,
            seq_from_str(
                "cgcagccatgacccagtcacgtagcgatagcggagtgtatactggcttaactatgcggcatcagagcagattgtactgagagtgcaccatatgcggtgtgaaataccgcacagatgcgtaaggagaaaataccgcatcagg",
            ),
        ),
        FeatureMapItem::new("BOM", Generic, seq_from_str("tttgtttaactttaagaaggaga")),
        FeatureMapItem::new(
            "BOM",
            Generic,
            seq_from_str(
                "cctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcatctggtgcactctcagtacaatctgctctgatgccgcatagttaagccagtatacactccgctatcgctacgtgactgggtcatggctgcg",
            ),
        ),
    ];

    let mut result = Vec::new();

    for item in items {
        let (matches_fwd, matches_rev) = match_subseq(&item.seq, seq);

        for (i, matches) in [matches_fwd, matches_rev].into_iter().enumerate() {
            // todo: Consider short descriptive notes for each.
            for range in matches {
                let dir = if i == 0 {
                    FeatureDirection::Forward
                } else {
                    FeatureDirection::Reverse
                };
                let direction = match item.feature_type {
                    CodingRegion | Ori | Promoter | AntibioticResistance => dir,
                    _ => FeatureDirection::None,
                };
                result.push(Feature {
                    range,
                    feature_type: item.feature_type,
                    label: item.name.clone(),
                    direction,
                    ..Default::default()
                });
            }
        }
    }

    result
}

pub fn find_features(seq: &[Nucleotide]) -> Vec<Feature> {
    let mut result = Vec::new();

    result.append(&mut find_his_tags(seq));
    result.append(&mut find_misc(seq));

    result
}
