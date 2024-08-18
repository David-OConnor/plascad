use std::{
    cmp::{min, PartialOrd},
    collections::HashSet,
    fmt, io,
    io::ErrorKind,
    ops::RangeInclusive,
};

use bincode::{Decode, Encode};
use eframe::egui::{pos2, Pos2};

use crate::{
    gui::seq_view::{NT_WIDTH_PX, SEQ_ROW_SPACING_PX, TEXT_X_START, TEXT_Y_START},
    sequence::{seq_complement, Nucleotide},
    Color, State,
};

/// A replacement for std::RangeInclusive, but copy type, and directly-accessible (mutable) fields.
/// An official replacement is eventually coming, but not for a while likely.
#[derive(Clone, Copy, Debug, PartialEq, Encode, Decode)]
pub struct RangeIncl {
    pub start: usize,
    pub end: usize,
}

impl RangeIncl {
    pub fn new(start: usize, end: usize) -> Self {
        Self { start, end }
    }
    pub fn contains(&self, val: usize) -> bool {
        val >= self.start && val <= self.end
    }

    /// This function handles both the +1 nature of our range indexing,
    /// and error handling for out of bounds. (The latter of which panics at runtime)
    pub fn index_seq<'a, T>(&self, seq: &'a [T]) -> Option<&'a [T]> {
        if self.start < 1 || self.end + 1 > seq.len() {
            return None;
        }

        Some(&seq[self.start - 1..=self.end - 1])
    }

    pub fn len(&self) -> usize {
        self.end - self.start + 1
    }
}

impl fmt::Display for RangeIncl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}..{}", self.start, self.end)
    }
}

/// Utility function to linearly map an input value to an output
pub fn map_linear(val: f32, range_in: (f32, f32), range_out: (f32, f32)) -> f32 {
    // todo: You may be able to optimize calls to this by having the ranges pre-store
    // todo the total range vals.
    let portion = (val - range_in.0) / (range_in.1 - range_in.0);

    portion * (range_out.1 - range_out.0) + range_out.0
}

/// We use this for dividing a nucleotied sequence into rows, for display in a canvas UI.
/// Each range is a nucleotide index, using our 1-based system.
pub fn get_row_ranges(len: usize, chars_per_row: usize) -> Vec<RangeIncl> {
    let mut result = Vec::new();

    // todo: Round etc instead of adding 1?
    let num_rows = len / chars_per_row + 1; // todo: +/-1 etc?

    for row_i in 0..num_rows {
        let end = row_i * chars_per_row + chars_per_row;
        // let end = min(end, len + 1);
        let end = min(end, len);

        // Adjustments are per our 1-based indexing system.
        result.push(RangeIncl::new(row_i * chars_per_row + 1, end));
    }

    result
}

// todo: We currently don't use this as a standalone fn; wrap back into `seq_i_to_pixel` a/r.
/// Maps sequence index, as displayed on a manually-wrapped UI display, to row and column indices.
fn seq_i_to_col_row(seq_i: usize, row_ranges: &[RangeIncl]) -> (usize, usize) {
    let mut row = 0;
    let mut row_range = RangeIncl::new(0, 10);

    for (row_, range) in row_ranges.iter().enumerate() {
        if range.contains(seq_i) {
            row = row_;
            row_range = range.clone();
            break;
        }
    }

    let col = seq_i - row_range.start;

    (col, row)
}

/// Maps sequence index, as displayed on a manually-wrapped UI display, to the relative pixel.
pub fn seq_i_to_pixel(seq_i: usize, row_ranges: &[RangeIncl]) -> Pos2 {
    let (col, row) = seq_i_to_col_row(seq_i, row_ranges);
    // This adjustment is used for placing the cursor at position 0; prior to the first nucleotide.
    let col = if seq_i == 0 { -1. } else { col as f32 };

    pos2(
        TEXT_X_START + col * NT_WIDTH_PX,
        TEXT_Y_START + row as f32 * SEQ_ROW_SPACING_PX,
    )
}

pub fn pixel_to_seq_i(pixel: Pos2, row_ranges: &[RangeIncl]) -> Option<usize> {
    let row = ((pixel.y - TEXT_Y_START) / SEQ_ROW_SPACING_PX) as usize;
    let col = ((pixel.x - TEXT_X_START) / NT_WIDTH_PX) as usize;

    // todo: Index vice loop?
    for (row_, range) in row_ranges.iter().enumerate() {
        if row_ == row {
            return Some(range.start + col);
        }
    }

    None
}

// todo; Move to Util A/R
pub fn remove_duplicates<T: Eq + std::hash::Hash>(vec: Vec<T>) -> Vec<T> {
    let set: HashSet<_> = vec.into_iter().collect();
    let vec: Vec<_> = set.into_iter().collect();
    vec
}

/// Given an index range of a feature, return sequence ranges for each row the feature occupies, that
/// contain the sequence. This, after converting to pixels, corresponds to how we draw features and primers.
/// This is used to draw overlays over the sequence that line up with a given index range.
/// todo: This should be RangeInclusive.
pub fn get_feature_ranges(feature_rng: &RangeIncl, all_ranges: &[RangeIncl]) -> Vec<RangeIncl> {
    let mut result = Vec::new();

    if feature_rng.end < feature_rng.start || feature_rng.end == 0 {
        eprintln!(
            "Error with feature ranges; start after end. Start: {}, end: {}",
            feature_rng.start, feature_rng.end
        );
        return result;
    }

    for range in all_ranges {
        if range.end < range.start || range.end == 0 {
            // eprintln!(
            //     "Error with ranges; start after end. Start: {}, end: {}",
            //     range.start,
            //     range.end
            // );
            return result;
        }

        if range.contains(feature_rng.start) {
            // Contains start only, or start and end.
            let end = min(range.end, feature_rng.end);

            // result.push(feature_rng.start..end - 1);
            result.push(RangeIncl::new(feature_rng.start, end));
        } else if range.contains(feature_rng.end) {
            // Contains end only.
            result.push(RangeIncl::new(range.start, feature_rng.end));
        } else if feature_rng.start < range.start && feature_rng.end > range.end {
            // Include this entire row
            result.push(RangeIncl::new(range.start, range.end));
        }
        // If none of the above, this row doesn't contain any of the sequence of interest.
    }

    result
}

pub fn color_from_hex(hex: &str) -> io::Result<Color> {
    let hex = &hex[1..]; // Remove the leading #

    if hex.len() < 6 {
        return Err(io::Error::new(ErrorKind::InvalidData, "Invalid color"));
    }

    let r = u8::from_str_radix(&hex[0..2], 16)
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Invalid color radix"))?;
    let g = u8::from_str_radix(&hex[2..4], 16)
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Invalid color radix"))?;
    let b = u8::from_str_radix(&hex[4..6], 16)
        .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Invalid color radix"))?;

    Ok((r, g, b))
}

pub fn color_to_hex(color: Color) -> String {
    format!("#{:x}{:x}{:x}", color.0, color.1, color.2)
}

/// Change the origin. This involves updating the sequence, and all features.
pub fn change_origin(state: &mut State) {
    let origin = &state.ui.new_origin;
    // Note the 1-based indexing logic we use.
    if *origin < 1 || *origin > state.generic.seq.len() {
        return;
    }

    state.generic.seq.rotate_left(origin - 1);

    let seq_len = state.generic.seq.len();
    for feature in &mut state.generic.features {
        // Convert to i32 to prevent an underflow on crash if we wrap. Use `rem_euclid`,
        // as Rust has unexpected behavior when using modulus on negatives.
        feature.range = RangeIncl::new(
            (feature.range.start as i32 + 1 - *origin as i32).rem_euclid(seq_len as i32) as usize,
            (feature.range.end as i32 + 1 - *origin as i32).rem_euclid(seq_len as i32) as usize,
        )
    }

    // todo: What else to update?
    state.sync_seq_related(None);
}

/// Find indexes where a subsequence matches a larger one, in both directions. Can be used to match primers,
/// known sequences etc.
/// Note: the reverse indices are reversed.
///         // todo: Partial matches as well.
pub fn match_subseq(subseq: &[Nucleotide], seq: &[Nucleotide]) -> (Vec<RangeIncl>, Vec<RangeIncl>) {
    let mut result = (Vec::new(), Vec::new()); // Forward, reverse

    let seq_len = seq.len();
    let subseq_len = subseq.len();
    let complement = seq_complement(seq);

    for seq_start in 0..seq_len {
        // Note: This approach handles sequence wraps, eg [circular] plasmids.
        let seq_iter = seq.iter().cycle().skip(seq_start).take(subseq_len);
        if subseq.iter().eq(seq_iter) {
            let seq_end = (seq_start + subseq_len) % seq_len;
            result.0.push(RangeIncl::new(seq_start, seq_end));
        }
    }

    for seq_start in 0..seq_len {
        let seq_iter = complement.iter().cycle().skip(seq_start).take(subseq_len);
        if subseq.iter().eq(seq_iter) {
            let seq_end = (seq_start + subseq_len) % seq_len;
            result.1.push(RangeIncl::new(seq_start, seq_end));
        }
    }

    result
}

// use std::env::current_exe;
// use winreg::enums::HKEY_CURRENT_USER;
// use winreg::RegKey;
//
// /// Associate file extensions in Windows. Requires the program to be run as an administrator.
// fn associate_windows_file_extension() -> io::Result<()> {
//     let hkcu = RegKey::predef(HKEY_CURRENT_USER);
//
//     let extension_key = hkcu.create_subkey(r"Software\Classes\.pcad")?;
//     extension_key.set_value("", &"pcadfile")?;
//
//     let file_type_key = hkcu.create_subkey(r"Software\Classes\pcadfile")?;
//     file_type_key.set_value("", &"My PCAD File")?;
//
//     let exe_path = current_exe()?.to_str().unwrap().to_string();
//
//     let command_key = hkcu.create_subkey(r"Software\Classes\pcadfile\shell\open\command")?;
//     command_key.set_value("", &format!("\"{}\" \"%1\"", exe_path))?;
//
//     Ok(())
// }
