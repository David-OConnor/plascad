use crate::{
    gui::navigation::{Page, PageSeq},
    primer::make_cloning_primers,
    sequence::{Feature, FeatureDirection, FeatureType, Seq},
    util::RangeIncl,
    Selection, State,
};

#[derive(Default)]
pub struct CloningInsertData {
    /// We use this list to store a list of features to clone an insert from, loaded from a file.
    pub features_loaded: Vec<Feature>,
    /// `cloning_ins_features_loaded` indexes reference this sequence.
    pub seq_loaded: Seq,
    pub feature_selected: Option<usize>,
    pub seq_insert: Seq,
    pub seq_input: String,
}

/// Create a new tab containing of the cloning product.
pub fn make_product_tab(state: &mut State) {
    // Note: This segment is almost a duplicate of `State::add_tab`, but retaining the generic data.
    state.generic.push(state.generic[state.active].clone());
    state
        .ion_concentrations
        .push(state.ion_concentrations[state.active].clone());
    state.path_loaded.push(None);
    state.portions.push(Default::default());

    state.active = state.generic.len() - 1;
    // state.add_tab();

    // Make sure to create cloning primers before performing the insert, or the result will be wrong.
    make_cloning_primers(state);

    // todo: Unecessary clone
    let insert = state.ui.cloning_insert.seq_insert.clone();
    state.insert_nucleotides(&insert, state.cloning_insert_loc);

    let label = match state.ui.cloning_insert.feature_selected {
        Some(i) => state.ui.cloning_insert.features_loaded[i].label.clone(),
        None => "Cloning insert".to_owned(),
    };

    // todo: Eventually, we'll likely be pulling in sequences already associated with a feature;
    // todo: Use the already existing data instead.
    state.generic[state.active].features.push(Feature {
        range: RangeIncl::new(
            state.cloning_insert_loc + 1,
            state.cloning_insert_loc + state.ui.cloning_insert.seq_insert.len(),
        ),
        label,
        feature_type: FeatureType::CodingRegion,
        direction: FeatureDirection::Forward,
        ..Default::default()
    });

    state.generic[state.active].metadata.plasmid_name = "Cloning product".to_owned();

    state.ui.page = Page::Map;
    state.ui.page_seq = PageSeq::View;
    state.ui.selected_item = Selection::Feature(state.generic[state.active].features.len() - 1);
}
