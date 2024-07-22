mod pcr;
mod portions;
pub mod primer;
pub mod seq_view; // pub for a few consts

use bincode::{Decode, Encode};
use eframe::{
    egui,
    egui::{Color32, Context, Key, RichText, ScrollArea, Ui},
};

use crate::{util::save, State};

pub const WINDOW_WIDTH: f32 = 1300.;
pub const WINDOW_HEIGHT: f32 = 800.;

pub const WINDOW_TITLE: &str = "Plasmid tools";

pub const ROW_SPACING: f32 = 22.;
pub const COL_SPACING: f32 = 30.;

#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum Page {
    /// Primer design and QC, including for cloning
    Primers,
    /// Determine optimal PCR parameters
    Pcr,
    Portions,
    // Sequence,
    // Enzymes,
    // Features,
}

impl Default for Page {
    fn default() -> Self {
        Self::Primers
    }
}

impl Page {
    pub fn to_str(self) -> String {
        match self {
            Self::Primers => "Primers",
            Self::Pcr => "PCR",
            Self::Portions => "Mixing portions",
        }
        .to_owned()
    }
}

#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum PagePrimerCreation {
    /// Eg amplifying a section of a single sequence.
    Amplification,
    SlicFc,
}

impl Default for PagePrimerCreation {
    fn default() -> Self {
        Self::SlicFc
    }
}

impl PagePrimerCreation {
    pub fn to_str(self) -> String {
        match self {
            Self::Amplification => "Amplification",
            Self::SlicFc => "SLIC and FastCloning",
        }
        .to_owned()
    }
}

#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum PageSeq {
    Edit,
    View,
}

impl Default for PageSeq {
    fn default() -> Self {
        Self::View
    }
}

impl PageSeq {
    pub fn to_str(self) -> String {
        match self {
            Self::Edit => "Edit squences",
            Self::View => "View sequences",
        }
        .to_owned()
    }
}

// // todo: Move A/R
// /// Used to determine which side of a primer we can extend or remove from in order to optimize it.
// #[derive(Clone, Copy, PartialEq)]
// enum TunableEnd {
//     None,
//     Left,
//     Right,
//     Both,
// }

// impl Default for TunableEnd {
//     fn default() -> Self {
//         Self::None
//     }
// }

fn page_button(page_state: &mut Page, page: Page, ui: &mut Ui) {
    let color = if *page_state == page {
        Color32::GREEN
    } else {
        Color32::WHITE
    };

    if ui
        .button(RichText::new(page.to_str()).color(color))
        .clicked()
    {
        *page_state = page;
    }

    ui.add_space(COL_SPACING);
}

fn page_selector(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        page_button(&mut state.ui.page, Page::Primers, ui);
        page_button(&mut state.ui.page, Page::Pcr, ui);
        page_button(&mut state.ui.page, Page::Portions, ui);
    });
}

fn page_primers_button(page_state: &mut PagePrimerCreation, page: PagePrimerCreation, ui: &mut Ui) {
    let color = if *page_state == page {
        Color32::GREEN
    } else {
        Color32::WHITE
    };

    if ui
        .button(RichText::new(page.to_str()).color(color))
        .clicked()
    {
        *page_state = page;
    }

    ui.add_space(COL_SPACING / 2.);
}

// todo: Use to_string and partialEq traits instead of duplicating the other page.
fn page_seq_button(page_state: &mut PageSeq, page: PageSeq, ui: &mut Ui) {
    let color = if *page_state == page {
        Color32::GREEN
    } else {
        Color32::WHITE
    };

    if ui
        .button(RichText::new(page.to_str()).color(color))
        .clicked()
    {
        *page_state = page;
    }

    ui.add_space(COL_SPACING / 2.);
}

pub fn page_primers_selector(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        page_primers_button(
            &mut state.ui.page_primer_creation,
            PagePrimerCreation::Amplification,
            ui,
        );
        page_primers_button(
            &mut state.ui.page_primer_creation,
            PagePrimerCreation::SlicFc,
            ui,
        );
    });
}

pub fn page_seq_selector(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        page_seq_button(&mut state.ui.page_seq, PageSeq::Edit, ui);
        page_seq_button(&mut state.ui.page_seq, PageSeq::View, ui);
    });
}

pub fn draw(state: &mut State, ctx: &Context) {
    ctx.input(|ip| {
        if ip.key_pressed(Key::A) && ip.modifiers.ctrl {
            state.primer_data.push(Default::default());
        }

        if ip.key_pressed(Key::S) && ip.modifiers.ctrl {
            if let Err(e) = save("plasmid_tools.save", state) {
                println!("Error saving: {e}");
            }
        }
    });

    egui::CentralPanel::default().show(ctx, |ui| {
        let mut visuals = ctx.style().visuals.clone();
        // visuals.override_text_color = Some(Color32::from_rgb(255, 0, 0));
        visuals.override_text_color = Some(Color32::LIGHT_GRAY);
        ctx.set_visuals(visuals);

        page_selector(state, ui);

        ui.add_space(ROW_SPACING);

        ScrollArea::vertical().show(ui, |ui| match state.ui.page {
            Page::Primers => primer::primer_page(state, ui),
            Page::Pcr => pcr::pcr_page(state, ui),
            Page::Portions => portions::portions_page(state, ui),
        });
    });
}
