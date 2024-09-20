use crate::{
    file_io::GenericData,
    gui::navigation::{Page, PageSeq},
    primer::make_cloning_primers,
    sequence::{seq_to_str, Feature, FeatureDirection, FeatureType, Seq},
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

/// Given a set of features and the sequence their ranges map to, set up our
/// insert sequences.
pub fn setup_insert_seqs(state: &mut State, features: Vec<Feature>, seq: Seq) {
    // todo: Unecessary cloning if loading from file.
    state.ui.cloning_insert.features_loaded = features;
    state.ui.cloning_insert.seq_loaded = seq;

    // Choose the initial insert as the CDS or gene with the largest len.
    let mut best = None;
    let mut best_len = 0;
    for (i, feature) in state.ui.cloning_insert.features_loaded.iter().enumerate() {
        let len = feature.len(state.generic[state.active].seq.len());
        if (feature.feature_type == FeatureType::CodingRegion
            || feature.feature_type == FeatureType::Gene)
            && len > best_len
        {
            best_len = len;
            best = Some(i);
        }
    }

    if let Some(feat_i) = best {
        let feature = &state.ui.cloning_insert.features_loaded[feat_i];

        if let Some(seq_this_ft) = feature.range.index_seq(&state.ui.cloning_insert.seq_loaded) {
            state.ui.cloning_insert.feature_selected = best;
            state.ui.cloning_insert.seq_insert = seq_this_ft.to_owned();
            state.ui.cloning_insert.seq_input = seq_to_str(seq_this_ft);
        }
    }
}

/// Create a new tab containing of the cloning product.
/// Optionally allow passing a new set of generic data to use, eg a backbone. If not present,
/// the current tab's will be used.
pub fn make_product_tab(state: &mut State, generic: Option<GenericData>) {
    // Note: This segment is almost a duplicate of `State::add_tab`, but retaining the generic data.

    let generic = match generic {
        Some(gen) => gen,
        None => state.generic[state.active].clone(),
    };

    state.generic.push(generic);
    state
        .ion_concentrations
        .push(state.ion_concentrations[state.active].clone());
    state.path_loaded.push(None);
    state.portions.push(Default::default());
    state.volatile.push(Default::default());

    state.active = state.generic.len() - 1;

    // Make sure to create cloning primers before performing the insert, or the result will be wrong.
    make_cloning_primers(state);

    // todo: Unecessary clone? Due to borrow rules.
    let insert = &state.ui.cloning_insert.seq_insert.clone();
    state.insert_nucleotides(insert, state.cloning_insert_loc);

    let label = match state.ui.cloning_insert.feature_selected {
        Some(i) => state.ui.cloning_insert.features_loaded[i].label.clone(),
        None => "Cloning insert".to_owned(),
    };

    // todo: Eventually, we'll likely be pulling in sequences already associated with a feature;
    // todo: Use the already existing data instead.
    state.generic[state.active].features.push(Feature {
        range: RangeIncl::new(
            state.cloning_insert_loc,
            state.cloning_insert_loc + state.ui.cloning_insert.seq_insert.len() - 1,
        ),
        label,
        feature_type: FeatureType::CodingRegion,
        direction: FeatureDirection::Forward,
        ..Default::default()
    });

    "Cloning product".clone_into(&mut state.generic[state.active].metadata.plasmid_name);

    state.ui.page = Page::Map;
    state.ui.page_seq = PageSeq::View;
    state.ui.selected_item = Selection::Feature(state.generic[state.active].features.len() - 1);
}
