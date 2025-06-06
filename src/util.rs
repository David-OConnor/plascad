use std::{cmp::min, collections::HashSet, fmt, io, io::ErrorKind, path::Path};

use bincode::{Decode, Encode};
use eframe::egui::{Pos2, pos2};
use na_seq::{
    Nucleotide,
    ligation::{filter_multiple_seqs, filter_unique_cutters, find_common_res},
    restriction_enzyme::RestrictionEnzyme,
    seq_complement,
};

use crate::{
    Color, ReUi,
    file_io::save::QUICKSAVE_FILE,
    gui::{
        WINDOW_TITLE,
        sequence::seq_view::{NT_WIDTH_PX, SEQ_ROW_SPACING_PX, TEXT_X_START, TEXT_Y_START},
    },
    misc_types::Feature,
    state::{State, StateVolatile},
};

const FEATURE_ANNOTATION_MATCH_THRESH: f32 = 0.95;

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
    pub fn index_seq<'a, T: Clone>(&self, seq: &'a [T]) -> Option<&'a [T]> {
        // todo: If end > len, wrap, and handle as circular.
        // if self.start < 1 || self.end > seq.len() {
        if self.start < 1 || self.end > seq.len() || self.start > self.end {
            return None;
        }

        // todo: This unfortunately doesn't work due to borrow rules.
        // We assume a circular sequence, for now. todo: Add as a parameter.
        // if self.start > self.end {
        //     let part1 = &seq[self.start - 1..];
        //     let part2 = &seq[..self.end - 1];
        //     Some(&[part1, part2].concat())
        // } else {
        //     Some(&seq[self.start - 1..=self.end - 1])
        // }

        Some(&seq[self.start - 1..=self.end - 1])
    }

    /// A wrap-based extension to `index_seq`.
    pub fn index_seq_wrap<T: Clone>(&self, seq: &[T]) -> Option<Vec<T>> {
        if self.start < 1 || self.end > seq.len() || self.start <= self.end {
            return None;
        }

        let part1 = seq[self.start - 1..].to_vec();
        // Todo: Not sure why we're including the end here without -1. Using for PCR to get it to come out right.
        let part2 = seq[..self.end].to_vec();
        let mut result = part1;
        result.extend(part2);
        Some(result)
    }

    pub fn len(&self) -> usize {
        if self.end < self.start {
            return 0;
        }
        self.end - self.start + 1
    }
}

impl fmt::Display for RangeIncl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}..{} {}bp", self.start, self.end, self.len())
    }
}

/// We use this for dividing a nucleotied sequence into rows, for display in a canvas UI.
/// Each range is a nucleotide index, using our 1-based system.
pub fn get_row_ranges(len: usize, chars_per_row: usize) -> Vec<RangeIncl> {
    // todo: Round etc instead of adding 1?
    let num_rows = len / chars_per_row + 1; // todo: +/-1 etc?

    let mut result = Vec::with_capacity(num_rows);

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
            row_range = *range;
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

//
/// For the sequence editor. Given an index range of a feature, return sequence ranges for each row the feature occupies, that
/// contain the sequence. This, after converting to pixels, corresponds to how we draw features and primers.
/// This is used to draw overlays over the sequence that line up with a given index range.
pub fn get_feature_ranges(
    feature_rng: &RangeIncl,
    all_ranges: &[RangeIncl],
    seq_len: usize,
) -> Vec<RangeIncl> {
    let mut result = Vec::new();

    // If the feature range wraps the origin, divide it into two ranges, and match both.
    let feature_ranges = if feature_rng.end < feature_rng.start {
        vec![
            RangeIncl::new(1, feature_rng.end),
            RangeIncl::new(feature_rng.start, seq_len),
        ]
    } else {
        vec![*feature_rng]
    };

    for ft_rng in &feature_ranges {
        for range in all_ranges {
            if range.end < range.start || range.end == 0 {
                eprintln!("Error with all ranges when finding a feature range: {range}");
                return result;
            }

            if range.contains(ft_rng.start) {
                // Contains start only, or start and end.
                let end = min(range.end, ft_rng.end);

                // result.push(ft_rng.start..end - 1);
                result.push(RangeIncl::new(ft_rng.start, end));
            } else if range.contains(ft_rng.end) {
                // Contains end only.
                result.push(RangeIncl::new(range.start, ft_rng.end));
            } else if ft_rng.start < range.start && ft_rng.end > range.end {
                // Include this entire row
                result.push(RangeIncl::new(range.start, range.end));
            }
            // If none of the above, this row doesn't contain any of the sequence of interest.
        }
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
    if *origin < 1 || *origin > state.get_seq().len() {
        return;
    }

    state.generic[state.active].seq.rotate_left(origin - 1);

    let seq_len = state.get_seq().len();
    for feature in &mut state.generic[state.active].features {
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
/// known sequences etc. Range indicies are relative to the forward direction.
/// todo: Partial matches as well.
pub fn match_subseq(subseq: &[Nucleotide], seq: &[Nucleotide]) -> (Vec<RangeIncl>, Vec<RangeIncl>) {
    let mut result = (Vec::new(), Vec::new()); // Forward, reverse

    let seq_len = seq.len();
    let subseq_len = subseq.len();
    let complement = seq_complement(seq);

    for seq_start in 0..seq_len {
        // Note: This approach handles sequence wraps, eg [circular] plasmids.
        let seq_iter = seq.iter().cycle().skip(seq_start).take(subseq_len);

        // let seq_iter: Vec<Nucleotide> = seq.iter().cycle().skip(seq_start).take(subseq_len).map(|nt| *nt).collect();

        // let similarity = seq_similarity(&seq_iter, subseq);
        // if similarity > FEATURE_ANNOTATION_MATCH_THRESH {
        //     let seq_end = (seq_start + subseq_len) % seq_len;
        //     result.0.push(RangeIncl::new(seq_start, seq_end));
        // }

        // println!("Similarity: {:?}", similarity);

        if subseq.iter().eq(seq_iter) {
            let seq_end = (seq_start + subseq_len) % seq_len;
            result.0.push(RangeIncl::new(seq_start + 1, seq_end));
        }
    }

    for seq_start in 0..seq_len {
        let seq_iter = complement.iter().cycle().skip(seq_start).take(subseq_len);
        // let seq_iter: Vec<Nucleotide> = complement.iter().cycle().skip(seq_start).take(subseq_len).map(|nt| *nt).collect();

        // let similarity = seq_similarity(&seq_iter, subseq);
        // if similarity > FEATURE_ANNOTATION_MATCH_THRESH {
        //     let seq_end = (seq_start + subseq_len) % seq_len;
        //     result.0.push(RangeIncl::new(seq_start, seq_end));
        // }
        //

        if subseq.iter().eq(seq_iter) {
            let seq_end = (seq_start + subseq_len) % seq_len;
            if seq_end < 1 {
                continue;
            }

            result
                .1
                .push(RangeIncl::new(seq_len - seq_end + 1, seq_len - seq_start));
        }
    }

    result
}

/// Find the similarity between the two sequences, on a scale of 0 to 1. Assumes same direction.
/// Note: This does not have a good way of handling length mismatches.
pub fn _seq_similarity(seq_a: &[Nucleotide], seq_b: &[Nucleotide]) -> f32 {
    // Using seq_a's len has consequences
    let mut matches = 0;
    for i in 0..seq_a.len() {
        if seq_a[i] == seq_b[i] {
            matches += 1;
        }
    }

    matches as f32 / seq_a.len() as f32
}

/// Merge a new set into an existing one. Don't add duplicates.
pub fn merge_feature_sets(existing: &mut Vec<Feature>, new: &[Feature]) {
    for feature_new in new {
        let mut exists = false;
        for feature_existing in &mut *existing {
            if feature_new.range == feature_existing.range {
                exists = true;
            }
        }
        if !exists {
            existing.push(feature_new.clone());
        }
    }
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

/// Get the title to be displayed in the windows tilebar.
pub fn get_window_title(path: &Path) -> String {
    let filename = path
        .file_name()
        .and_then(|name| name.to_str())
        .map(|name_str| name_str.to_string())
        .unwrap();

    if filename == QUICKSAVE_FILE {
        WINDOW_TITLE.to_owned()
    } else {
        filename
    }
}

/// We filter for restriction enzymes based on preferences set. We do this in several stages.
/// Note that this function includes filter characteristics that inolve matches across
/// multiple opened tabs (sequences).
pub fn filter_res<'a>(
    data: &ReUi,
    volatile: &[StateVolatile],
    lib: &'a [RestrictionEnzyme],
) -> Vec<&'a RestrictionEnzyme> {
    let mut re_match_set = Vec::new(); // By tab
    for active in &data.tabs_selected {
        re_match_set.push(&volatile[*active].restriction_enzyme_matches);
    }

    let mut result = find_common_res(&re_match_set, lib, data.sticky_ends_only);

    if data.multiple_seqs {
        filter_multiple_seqs(&mut result, &re_match_set, lib);
    }

    // If `unique_cutters_only` is selected, each RE must be unique in each sequence. If it's two or more times in any given sequence, don't
    // add it.
    if data.unique_cutters_only {
        filter_unique_cutters(&mut result, &re_match_set, lib);
    }

    result
}
