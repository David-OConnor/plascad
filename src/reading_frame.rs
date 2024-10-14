use std::fmt::Display;

use bincode::{Decode, Encode};
use na_seq::{
    seq_complement, Nucleotide,
    Nucleotide::{A, G, T},
    Seq,
};

use crate::util::RangeIncl;

const START_CODON: [Nucleotide; 3] = [A, T, G];
pub const STOP_CODONS: [[Nucleotide; 3]; 3] = [[T, A, A], [T, A, G], [T, G, A]];

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

    pub fn is_reverse(&self) -> bool {
        // todo: Enum
        !matches!(self, Self::Fwd0 | Self::Fwd1 | Self::Fwd2)
    }

    /// Get a seqeuence of the full sequence in the appropriate direction, and with the appropriate offset.
    /// The result will always start in frame. It will be trimmed if offset > 0.
    /// todo: THis makes a clone. Can we instead do a slice?
    pub fn arrange_seq(&self, seq: &[Nucleotide]) -> Seq {
        let offset = self.offset();

        if self.is_reverse() {
            seq_complement(seq)[offset..].to_vec()
        } else {
            seq[offset..].to_vec()
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

#[derive(Debug, Clone, Encode, Decode)]
pub struct ReadingFrameMatch {
    pub frame: ReadingFrame,
    /// Indices are respective to the non-complementary seq, for both forward and reverse reading frames.
    pub range: RangeIncl,
}

/// Find coding regions in a sequence, given a reading frame.
pub fn find_orf_matches(seq: &[Nucleotide], orf: ReadingFrame) -> Vec<ReadingFrameMatch> {
    let mut result = Vec::new();

    let offset = orf.offset();
    let seq_len_full = seq.len();

    if seq_len_full < 3 {
        return result;
    }

    let seq_ = orf.arrange_seq(&seq);
    let len = seq_.len();

    let mut frame_open = None; // Inner: Start index.

    for i_ in 0..len / 3 {
        let i = i_ * 3; // The actual sequence index.

        let nts = &seq_[i..i + 3];

        if frame_open.is_none() && nts == START_CODON {
            frame_open = Some(i);
        // } else if frame_open.is_some() && stop_codons.contains(nts.try_into().unwrap()) {
        } else if frame_open.is_some()
            && (STOP_CODONSit commit .contains(nts.try_into().unwrap()) || seq_len_full - i <= 3)
        {
            // If we reach the end of the sequence, consider it closed.
            // todo: Handle circular around the origin. Ie, don't auto-close in that case.
            // + 1 for our 1-based seq name convention.
            // This section's a bit hairy; worked by trial and error. Final indices are respective to
            // the non-complementary seq, for both forward and reverse reading frames.
            let range = if !orf.is_reverse() {
                // todo: We still have wonkiness here.
                // RangeIncl::new(frame_open.unwrap() + 1 + offset, i + 2 + offset)
                RangeIncl::new(frame_open.unwrap() + 1 + offset, i + 3 + offset)
            } else {
                RangeIncl::new(
                    seq_len_full - (i + 2 + offset),
                    seq_len_full - (frame_open.unwrap() + offset) - 1,
                )
            };

            result.push(ReadingFrameMatch { frame: orf, range });
            frame_open = None;
        }
    }

    result
}
