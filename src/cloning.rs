//! Used for a semi-automated cloning process that chooses a suitable backbone and restriction enzymes
//! or primers.
//!
//! todo: Make sure you're QCing and filtering for expression system.
//!

//! todo: Allow users to enter custom backbones.

use crate::{
    backbones::Backbone,
    file_io::GenericData,
    gui::navigation::{Page, PageSeq},
    ligation::{filter_multiple_seqs, filter_unique_cutters, find_common_res},
    primer::make_cloning_primers,
    restriction_enzyme::{find_re_matches, RestrictionEnzyme},
    sequence::{seq_to_str, Feature, FeatureDirection, FeatureType, Nucleotide, Seq},
    util::RangeIncl,
    Selection, State, StateVolatile,
};

/// Include this many nucleotides to the left, and right of each insert, when searching for RE sites.
/// note: 4-6 nucleotides may be an ideal buffer. Note that this should be conservatively long.
pub const RE_INSERT_BUFFER: usize = 22;

/// Attempt to insert the sequence this far downstream of the RBS. 5-10 is ideal.
pub const RBS_BUFFER: usize = 7;

pub const RBS_BUFFER_MIN: isize = 4;
pub const RBS_BUFFER_MAX: isize = 11;

#[derive(Clone, Copy, PartialEq)]
pub enum Status {
    Pass,
    Fail,
    NotApplicable,
}

impl Default for Status {
    fn default() -> Self {
        Self::Fail
    }
}

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

/// For a given insert and vector, find suitable  restriction enzymes for cloning.
/// Make sure that the insert sequence is properly buffered to allow for RE matching outside
/// of the coding region (etc)'s range, upstream of this.
/// todo: Fow now, we limit our options to unique-match, single-cutters.
pub fn find_re_candidates<'a>(
    backbone: &Backbone,
    seq_insert: &[Nucleotide],
    // re_match_sets: &[&[ReMatch]],
    // insert_tab: usize,
    // insert_range: RangeIncl,
    lib: &'a [RestrictionEnzyme],
    volatile: &[StateVolatile],
    // ) -> Vec<&'a RestrictionEnzyme> {
) -> Vec<RestrictionEnzyme> {
    // Note: The first part of this function is similar to how we filter REs for digest on the digest/ligation page.

    let matches_insert = find_re_matches(seq_insert, lib);
    let matches_backbone = find_re_matches(&backbone.seq, lib);

    let re_match_set = [&matches_insert, &matches_backbone];

    // Set up our initial REs: Ones that match both the backbone, and insert. (anywhere, for now)
    let mut result = find_common_res(&re_match_set, lib, true);

    // For now, we are always filtering by these, as well as setting sticky ends only.
    filter_multiple_seqs(&mut result, &re_match_set, lib);
    filter_unique_cutters(&mut result, &re_match_set, lib);

    // Filter for REs that are at an appropriate place on the insert.

    result.into_iter().map(|r| r.clone()).collect()
}

#[derive(Default)]
pub struct AutocloneStatus {
    pub rbs_dist: Status,
    pub downstream_of_promoter: Status,
    pub upstream_of_terminator: Status,
    // todo: Check for terminator?
    pub direction: Status,
    pub tag_frame: Status,
    // pub primer_quality: Status,
    // pub re_dist: Status,
}

impl AutocloneStatus {
    pub fn new(backbone: &Backbone, insert_loc: usize, insert_len: usize) -> Self {
        let rbs_dist = match backbone.rbs {
            Some(rbs) => {
                let dist = insert_loc as isize - rbs.end as isize;
                // todo:  Handle wraps.
                if dist >= RBS_BUFFER_MIN && dist <= RBS_BUFFER_MAX {
                    Status::Pass
                } else {
                    Status::Fail
                }
            }
            None => Status::NotApplicable,
        };

        // todo: Handle wraps.
        let downstream_of_promoter = match backbone.promoter {
            Some(p) => {
                if insert_loc > p.end {
                    Status::Pass
                } else {
                    Status::Fail
                }
            }
            None => Status::NotApplicable,
        };

        let upstream_of_terminator = match backbone.terminator {
            Some(p) => {
                if insert_loc > p.end {
                    Status::Pass
                } else {
                    Status::Fail
                }
            }
            None => Status::NotApplicable,
        };

        let direction = Status::Pass; // todo

        let tag_frame = match backbone.his_tag {
            Some(tag_range) => {
                let dist_from_start_codon = (tag_range.start - insert_loc) + insert_len;
                if dist_from_start_codon % 3 == 0 {
                    Status::Pass
                } else {
                    Status::Fail
                }
            }
            None => Status::NotApplicable,
        };

        Self {
            rbs_dist,
            downstream_of_promoter,
            upstream_of_terminator,
            direction,
            tag_frame,
        }
    }
}
