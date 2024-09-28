// todo: unused; delete when ready.

use eframe::{
    egui::{pos2, vec2, Frame, Pos2, Rect, RichText, Sense, Shape, Stroke, TextEdit, Ui},
    emath::RectTransform,
    epaint::Color32,
};

use crate::{
    cloning::make_product_tab,
    gui::{
        autocloning,
        lin_maps::{lin_map_zoomed, OFFSET},
        BACKGROUND_COLOR, COL_SPACING, LINEAR_MAP_HEIGHT, ROW_SPACING,
    },
    sequence::{seq_from_str, seq_to_str},
    util::map_linear,
    State,
};

/// Draw a mini sequence display in its own canvas. Zoomed in on the cloning insert location, and with
/// a line drawn on it.
pub fn seq_lin_disp_cloning(state: &State, ui: &mut Ui) {
    Frame::canvas(ui.style())
        .fill(BACKGROUND_COLOR)
        .show(ui, |ui| {
            let (response, _painter) = {
                let desired_size = vec2(ui.available_width(), LINEAR_MAP_HEIGHT);
                ui.allocate_painter(desired_size, Sense::click())
            };

            let to_screen = RectTransform::from_to(
                Rect::from_min_size(Pos2::ZERO, response.rect.size()),
                response.rect,
            );

            const NT_WIDTH: usize = 400;

            let mut shapes = lin_map_zoomed(
                &state,
                &to_screen,
                state.cloning_insert_loc,
                NT_WIDTH,
                state.active,
                ui,
            );

            // todo: The DRY sections below are related to drawing the vertical line at the insert point.

            let seq_len = state.get_seq().len();

            if seq_len == 0 {
                return; // Avoid divide-by-0
            }

            let index_left = (state.cloning_insert_loc as isize - (NT_WIDTH / 2) as isize)
                .rem_euclid(seq_len as isize) as usize; // Rust awk % on negative values.
            let index_right = (state.cloning_insert_loc + NT_WIDTH / 2) % seq_len;
            let pixel_left = OFFSET.x;
            let pixel_right = ui.available_width() - 2. * OFFSET.x;

            // todo: This segment is DRY from draw_linear_map. Standalone fn etcx. {}    let index_to_x = |mut i: usize| {
            let index_to_x = |mut i: usize| {
                // This handles the case when the zoomed-in view is near the top; the left index
                // will be near the end of the sequence, incorrectly calculating the portion-through in
                // the linear map.
                let right = if index_left > index_right {
                    if i < index_right {
                        // Ie, we are to the right of the origin.
                        i += seq_len
                    }

                    index_right + seq_len
                } else {
                    index_right
                };

                map_linear(
                    i as f32,
                    (index_left as f32, right as f32),
                    (pixel_left, pixel_right),
                )
            };

            // Draw the insertion site.
            let point_top = pos2(index_to_x(state.cloning_insert_loc), 4.);
            let point_bottom = pos2(index_to_x(state.cloning_insert_loc), 44.);

            shapes.push(Shape::line_segment(
                [to_screen * point_bottom, to_screen * point_top],
                Stroke::new(3., Color32::YELLOW),
            ));

            ui.painter().extend(shapes);
        });
}

pub fn seq_editor_slic(state: &mut State, ui: &mut Ui) {
    ui.heading("SLIC and FastCloning");

    // todo: Once you add this capability.
    ui.label("Clone a sequence into this one. Below, either paste the insert sequence, or select a \
    file (GenBank, SnapGene, FASTA, or PlasCAD) containing the insert sequence. Set the insert location: \
    This is the position in the vector (The currently open sequence) the insert will be placed after. Then click \"Clone\". This will \
    create and open a new sequence, and create optimized primers for both the insert and vector.");

    ui.add_space(ROW_SPACING);

    // ui.label(
    //     "A file dialog will open, prompting you to save a new file for the combined product. Your \
    // current (vector) file will be saved, and the new cloning product file will be opened.",
    // );

    // todo: Zoomed in on insert loc, and draw insert loc.
    seq_lin_disp_cloning(state, ui);
    ui.add_space(ROW_SPACING / 2.);

    ui.horizontal(|ui| {
        ui.label("Insert location: ");
        let mut entry = state.cloning_insert_loc.to_string();
        if ui
            .add(TextEdit::singleline(&mut entry).desired_width(40.))
            .changed()
        {
            state.cloning_insert_loc = entry.parse().unwrap_or(0);
        }

        ui.add_space(COL_SPACING);

        if state.ui.cloning_insert.seq_insert.len() > 6 {
            if ui
                .button(RichText::new("Clone").color(Color32::GOLD))
                .clicked()
            {
                // Save this vector; this file or quicksave instance will be turned into the cloning
                // product.
                make_product_tab(state, None);
            }
        }
    });

    ui.add_space(ROW_SPACING);
    autocloning::insert_file_section(state, ui);

    ui.add_space(ROW_SPACING);
    // Note: Unlike RE cloning, we don't want a buffer region, hence passing 0 here.
    autocloning::insert_selector(&mut state.ui.cloning_insert, 0, ui);

    ui.add_space(ROW_SPACING);
    ui.horizontal(|ui| {
        ui.heading("Insert:");
        ui.label(&format!(
            "len: {}",
            state.ui.cloning_insert.seq_insert.len()
        ));
    });

    ui.add_space(ROW_SPACING);
}
