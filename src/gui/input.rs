//! Code related to mouse and keyboard input handling.

use std::{mem, path::PathBuf};

use eframe::egui::{Event, InputState, Key, PointerButton, Ui};
use na_seq::{seq_from_str, Nucleotide};

use crate::{
    file_io::{
        save,
        save::{load_import, StateToSave, QUICKSAVE_FILE},
    },
    gui::{
        navigation::{Page, Tab},
        set_window_title,
    },
    state::State,
    util::RangeIncl,
    StateUi,
};

/// Handle hotkeys and clicks that affect all pages.
fn handle_global(state: &mut State, ip: &InputState) {
    if ip.key_pressed(Key::A) && ip.modifiers.ctrl && !state.ui.text_edit_active {
        if !state.get_seq().is_empty() {
            state.ui.text_selection = Some(RangeIncl::new(1, state.get_seq().len()))
        }
    }

    if ip.key_pressed(Key::S) && ip.modifiers.ctrl && !ip.modifiers.shift {
        save::save_current_file(state);
    }

    if ip.key_pressed(Key::N) && ip.modifiers.ctrl {
        state.add_tab();
        state.tabs_open.push(Tab {
            path: None,
            ab1: false,
        });
    }

    if ip.key_pressed(Key::F) && ip.modifiers.ctrl {
        state.ui.highlight_search_input = true;
        // Disable the cursor, so we don't insert nucleotides while searching!
        state.ui.text_cursor_i = None;
    }

    if ip.key_pressed(Key::S) && ip.modifiers.ctrl && ip.modifiers.shift {
        state.ui.file_dialogs.save.select_file();
    }

    if ip.key_pressed(Key::O) && ip.modifiers.ctrl {
        state.ui.file_dialogs.load.select_file();
    }

    state.ui.cursor_pos = ip.pointer.hover_pos().map(|pos| (pos.x, pos.y));

    if ip.pointer.button_clicked(PointerButton::Primary) {
        // todo: Not working for fixing our fast-sel off-by-one bug.
        // if let Some(i) = &state.ui.cursor_seq_i {
        //     state.ui.text_selection = Some(*i..=*i); // 1-based indexing. Second value is a placeholder.
        // }

        state.ui.click_pending_handle = true;
    }

    if ip.pointer.button_double_clicked(PointerButton::Primary) {
        state.ui.dblclick_pending_handle = true;
    }
}

/// Handle sequence selection on the sequence page, as when dragging the mouse.
fn handle_seq_selection(state_ui: &mut StateUi, dragging: bool) {
    if dragging {
        // A drag has started.
        if !state_ui.dragging {
            state_ui.dragging = true;
            // We are handling in the seq view after setting cursor seq i. Still glitchy.

            if let Some(i) = &state_ui.cursor_seq_i {
                state_ui.text_selection = Some(RangeIncl::new(*i, *i)); // The end value is a placeholder.
            }
        } else {
            // The drag is in progress; continually update the selection, for visual feedback.
            if let Some(i) = &state_ui.cursor_seq_i {
                if let Some(sel_range) = &mut state_ui.text_selection {
                    sel_range.end = *i;
                }
            }
        }
    } else {
        // A drag has ended.
        if state_ui.dragging {
            state_ui.dragging = false;

            // This logic allows for reverse-order selecting, ie dragging the cursor from the end to
            // the start of the region. Handle in the UI how to display this correctly while the drag is in process,
            // since this only resolves once it's complete.
            if let Some(sel_range) = &mut state_ui.text_selection {
                if sel_range.start > sel_range.end {
                    mem::swap(&mut sel_range.start, &mut sel_range.end);
                }
            }
        }
    }
}

/// Handles keyboard and mouse input not associated with a widget.
/// todo: MOve to a separate module if this becomes complex.
pub fn handle_input(state: &mut State, ui: &mut Ui) {
    let mut reset_window_title = false; // This setup avoids borrow errors.

    ui.ctx().input(|ip| {
        // Check for file drop
        if let Some(dropped_files) = ip.raw.dropped_files.first() {
            if let Some(path) = &dropped_files.path {
                if let Some(loaded) = load_import(path) {
                    state.load(&loaded);
                    reset_window_title = true;
                }
            }
        }

        handle_global(state, ip);

        if state.ui.text_edit_active && !reset_window_title {
            return;
        }

        // This event match is not specific to the seqe page
        for event in &ip.events {
            match event {
                Event::Cut => state.copy_seq(),
                Event::Copy => state.copy_seq(),
                Event::Paste(pasted_text) => {}
                _ => (),
            }
        }

        if let Page::Sequence = state.ui.page {
            // This is a bit awk; borrow errors.
            let mut move_cursor: Option<i32> = None;
            // todo: How can we control the rate?
            if state.ui.text_cursor_i.is_some() {
                if ip.key_pressed(Key::ArrowLeft) {
                    move_cursor = Some(-1);
                }
                if ip.key_pressed(Key::ArrowRight) {
                    move_cursor = Some(1);
                }
                if ip.key_pressed(Key::ArrowUp) {
                    move_cursor = Some(-(state.ui.nt_chars_per_row as i32));
                }
                if ip.key_pressed(Key::ArrowDown) {
                    move_cursor = Some(state.ui.nt_chars_per_row as i32);
                }
                // Escape key: Remove the text cursor.
                if ip.key_pressed(Key::Escape) {
                    move_cursor = None;
                    state.ui.text_cursor_i = None;
                    state.ui.text_selection = None;
                }
            }

            if ip.key_pressed(Key::Escape) {
                state.ui.text_selection = None;
            }

            // Insert nucleotides A/R.
            if let Some(mut i) = state.ui.text_cursor_i {
                if !state.ui.seq_edit_lock {
                    if i > state.get_seq().len() {
                        i = 0; // todo?? Having an overflow when backspacing near origin.
                    }

                    let i = i + 1; // Insert after this nucleotide; not before.
                                   // Don't allow accidental nt insertion when the user is entering into the search bar.

                    // Add NTs.
                    if ip.key_pressed(Key::A) && !ip.modifiers.ctrl {
                        state.insert_nucleotides(&[Nucleotide::A], i);
                    }
                    if ip.key_pressed(Key::T) {
                        state.insert_nucleotides(&[Nucleotide::T], i);
                    }
                    if ip.key_pressed(Key::C) && !ip.modifiers.ctrl {
                        state.insert_nucleotides(&[Nucleotide::C], i);
                    }
                    if ip.key_pressed(Key::G) {
                        state.insert_nucleotides(&[Nucleotide::G], i);
                    }
                    if ip.key_pressed(Key::Backspace) && i > 1 {
                        state.remove_nucleotides(RangeIncl::new(i - 2, i - 2));
                    }
                    if ip.key_pressed(Key::Delete) {
                        state.remove_nucleotides(RangeIncl::new(i - 1, i - 1));
                    }

                    // Paste nucleotides
                    for event in &ip.events {
                        match event {
                            Event::Cut => {
                                // state.remove_nucleotides();
                                // move_cursor = Some(pasted_text.len() as i32);
                            }
                            Event::Copy => {}
                            Event::Paste(pasted_text) => {
                                state.insert_nucleotides(&seq_from_str(pasted_text), i);
                                move_cursor = Some(pasted_text.len() as i32);
                            }
                            _ => (),
                        }
                    }
                }
            }

            if !state.ui.seq_edit_lock {
                if ip.key_pressed(Key::A) && !ip.modifiers.ctrl {
                    move_cursor = Some(1);
                }
                if ip.key_pressed(Key::T) {
                    move_cursor = Some(1);
                }
                if ip.key_pressed(Key::C) && !ip.modifiers.ctrl {
                    move_cursor = Some(1);
                }
                if ip.key_pressed(Key::G) {
                    move_cursor = Some(1);
                }

                if ip.key_pressed(Key::Backspace) {
                    move_cursor = Some(-1);
                }
            }

            if let Some(i) = &mut state.ui.text_cursor_i {
                if let Some(amt) = move_cursor {
                    let val = *i as i32 + amt;
                    // If val is 0, the cursor is before the first character; this is OK.
                    // If it's < 0, we will get an overflow when converting to usize.
                    if val >= 0 {
                        let new_posit = val as usize;
                        if new_posit <= state.generic[state.active].seq.len() {
                            *i = new_posit;
                        }
                    }
                }
            }

            handle_seq_selection(&mut state.ui, ip.pointer.is_decidedly_dragging());
        }
    });

    if reset_window_title {
        set_window_title(&state.tabs_open[state.active], ui);
    }
}
