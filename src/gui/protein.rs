use eframe::egui::{Color32, FontFamily, FontId, RichText, Ui};

use crate::{
    amino_acids::{AaIdent, AminoAcid},
    gui::{COL_SPACING, ROW_SPACING},
    sequence::FeatureType,
    State,
};

const COLOR_PROT_SEQ: Color32 = Color32::from_rgb(255, 100, 200);
const FONT_SIZE_SEQ: f32 = 14.;

// todo: Eval how cacheing and state is handled.

pub fn protein_page(state: &mut State, ui: &mut Ui) {
    ui.horizontal(|ui| {
        ui.heading("Proteins, from sequence coding regions");

        ui.add_space(COL_SPACING);
        ui.label("One letter ident:");
        let mut one_letter = state.ui.aa_ident_disp == AaIdent::OneLetter;
        if ui.checkbox(&mut one_letter, "").changed() {
            state.ui.aa_ident_disp = if one_letter {
                AaIdent::OneLetter
            } else {
                AaIdent::ThreeLetters
            };
        }
    });

    ui.add_space(ROW_SPACING);

    for (i, feature) in state.generic.features.iter().enumerate() {
        if feature.feature_type != FeatureType::CodingRegion {
            continue;
        }

        for ((j, om)) in &state.volatile.cr_orf_matches {
            if *j == i {
                ui.heading(RichText::new(&feature.label).color(Color32::LIGHT_BLUE));
                ui.label(format!("Reading frame: {}, Range: {}", om.frame, om.range));

                // todo: Consider switching to a canvas-based approach A/R. Text as an MVP.
                // todo: Breakout into fns, probably in src/protein.
                if let Some(seq_orf_match_dna) = om.range.index_seq(&state.generic.seq) {
                    // todo: DRy with feature_db_load.
                    let len = seq_orf_match_dna.len();

                    let mut aa_seq = Vec::new();

                    for i_ in 0..len / 3 {
                        let i = i_ * 3; // The actual sequence index.

                        let nts = &seq_orf_match_dna[i..i + 3];

                        let mut matched = false;
                        if let Some(aa) = AminoAcid::from_codons(nts.try_into().unwrap()) {
                            // todo: Handle unknown AAs; don't just silently ommit as we are currently doing.
                            aa_seq.push(aa);
                        }
                    }

                    let mut aa_text = String::new();
                    for aa in aa_seq {
                        let aa_str = match state.ui.aa_ident_disp {
                            AaIdent::ThreeLetters => &format!("{}  ", aa.ident_3_letter()),
                            AaIdent::OneLetter => &format!("{}  ", aa.ident_single_letter()),
                        };
                        aa_text.push_str(aa_str);
                    }
                    ui.label(
                        RichText::new(aa_text)
                            .color(COLOR_PROT_SEQ)
                            .font(FontId::new(FONT_SIZE_SEQ, FontFamily::Monospace)),
                    );
                }

                // let offset = om.frame.offset();
                // let seq_full = om.frame.arrange_seq(&state.generic.seq);
                // let len_full = seq_full.len();

                // todo: For now, constructing this dynamically. Cache it in state.volatile.
                // todo: Handle reverse reading frames.
                if let Some(seq_dna) = feature.range.index_seq(&state.generic.seq) {
                    // todo: DRy with feature_db_load.
                }

                ui.add_space(ROW_SPACING);
            }
        }
    }
}
