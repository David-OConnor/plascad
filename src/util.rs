use std::{cmp::min, collections::HashSet, num::ParseIntError, ops::Range};

use eframe::egui::{pos2, Pos2};

use crate::{
    gui::seq_view::{NT_WIDTH_PX, SEQ_ROW_SPACING_PX, TEXT_X_START, TEXT_Y_START},
    Color, State,
};
/// Utility function to linearly map an input value to an output
pub fn map_linear(val: f32, range_in: (f32, f32), range_out: (f32, f32)) -> f32 {
    // todo: You may be able to optimize calls to this by having the ranges pre-store
    // todo the total range vals.
    let portion = (val - range_in.0) / (range_in.1 - range_in.0);

    portion * (range_out.1 - range_out.0) + range_out.0
}

/// We use this for dividing a nucleotied sequence into rows, for display in a canvas UI.
pub fn get_row_ranges(len: usize, chars_per_row: usize) -> Vec<Range<usize>> {
    let mut result = Vec::new();

    // todo: Round etc instead of adding 1?
    let num_rows = len / chars_per_row + 1; // todo: +/-1 etc?

    for row_i in 0..num_rows {
        result.push(row_i * chars_per_row..row_i * chars_per_row + chars_per_row);
    }

    result
}

// todo: We currently don't use this as a standalone fn; wrap back into `seq_i_to_pixel` a/r.
/// Maps sequence index, as displayed on a manually-wrapped UI display, to row and column indices.
fn seq_i_to_col_row(seq_i: usize, row_ranges: &[Range<usize>]) -> (usize, usize) {
    let mut row = 0;
    let mut row_range = 0..10;

    for (row_, range) in row_ranges.iter().enumerate() {
        if range.contains(&seq_i) {
            row = row_;
            row_range = range.clone();
            break;
        }
    }

    let col = seq_i - row_range.start;

    (col, row)
}

/// Maps sequence index, as displayed on a manually-wrapped UI display, to the relative pixel.
pub fn seq_i_to_pixel(seq_i: usize, row_ranges: &[Range<usize>]) -> Pos2 {
    let (col, row) = seq_i_to_col_row(seq_i, row_ranges);

    pos2(
        TEXT_X_START + col as f32 * NT_WIDTH_PX,
        TEXT_Y_START + row as f32 * SEQ_ROW_SPACING_PX,
    )
}

pub fn pixel_to_seq_i(pixel: Pos2, row_ranges: &[Range<usize>]) -> Option<usize> {
    // todo: ROunding?
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

// todo: Abstract out/use in feature overlay
/// Given an index range of a feature, return sequence ranges for each row the feature occupies, that
/// contain the sequence. This, after converting to pixels, corresponds to how we draw features and primers.
pub fn get_feature_ranges(
    feature_rng: &Range<usize>,
    all_ranges: &[Range<usize>],
) -> Vec<Range<usize>> {
    let mut result = Vec::new();

    if feature_rng.end < feature_rng.start || feature_rng.end == 0 {
        // eprintln!("Error with feature ranges; start after end. Start: {}, end: {}", feature_rng.start, feature_rng.end);
        return result;
    }

    for range in all_ranges {
        if range.end < range.start || range.end == 0 {
            // eprintln!("Error with ranges; start after end. Start: {}, end: {}", range.start, range.end);
            return result;
        }

        if range.contains(&feature_rng.start) {
            // Contains start only, or start and end.
            let end = min(range.end, feature_rng.end);

            result.push(feature_rng.start..end - 1); // todo experimenting.
        } else if range.contains(&feature_rng.end) {
            // Contains end only.
            result.push(range.start..feature_rng.end); // todo experimenting.
        } else if feature_rng.start < range.start && feature_rng.end > range.end {
            // Include this entire row
            result.push(range.start..range.end - 1);
        }
        // If none of the above, this row doesn't contain any of the sequence of interest.
    }

    result
}

pub fn color_from_hex(hex: &str) -> Result<Color, ParseIntError> {
    // if hex.len() != 6 {
    //     return Err(std::num::ParseIntError::new()); // This line won't actually compile since there's no ParseIntError::new(), we'll handle the length check below.
    // }

    let hex = &hex[1..]; // Remove the leading #

    let r = u8::from_str_radix(&hex[0..2], 16)?;
    let g = u8::from_str_radix(&hex[2..4], 16)?;
    let b = u8::from_str_radix(&hex[4..6], 16)?;

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
        feature.index_range.0 =
            (feature.index_range.0 as i32 + 1 - *origin as i32).rem_euclid(seq_len as i32) as usize;
        feature.index_range.1 =
            (feature.index_range.1 as i32 + 1 - *origin as i32).rem_euclid(seq_len as i32) as usize;
    }

    // todo: What else to update?
    state.sync_seq_related(None);
}
