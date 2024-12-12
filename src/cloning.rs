//! Used for a semi-automated cloning process that chooses a suitable backbone and restriction enzymes
//! or primers.
//!
//! todo: Make sure you're QCing and filtering for expression system.
//!

//! todo: Allow users to enter custom backbones.

use na_seq::{
    insert_into_seq,
    ligation::{filter_multiple_seqs, filter_unique_cutters, find_common_res},
    restriction_enzyme::{find_re_matches, ReMatch, RestrictionEnzyme},
    seq_to_str_lower, AminoAcid, CodingResult, Nucleotide, Seq,
};

use crate::{
    backbones::Backbone,
    file_io::GenericData,
    gui::navigation::{Page, PageSeq},
    misc_types::{Feature, FeatureDirection, FeatureType},
    primer::{make_cloning_primers, Primer},
    state::State,
    util::RangeIncl,
    Selection,
};

/// Include this many nucleotides to the left, and right of each insert, when searching for RE sites.
/// note: 4-6 nucleotides may be an ideal buffer. Note that this should be conservatively long.
pub const RE_INSERT_BUFFER: usize = 22;

/// Attempt to insert the sequence this far downstream of the RBS. 5-10 is ideal.
pub const RBS_BUFFER: usize = 7;

pub const RBS_BUFFER_MIN: isize = 4;
pub const RBS_BUFFER_MAX: isize = 11;

pub struct CloningState {
    pub backbone_selected: BackboneSelected,
    /// Note: This is only used currently if using the opened file as the BB; otherwise
    /// we use a ref to the library.
    pub backbone: Option<Backbone>, // todo: Other options: Store a ref; use in bb_selected instead of an index.
    pub res_common: Vec<RestrictionEnzyme>,
    pub re_matches_vec_common: Vec<ReMatch>,
    pub re_matches_insert_common: Vec<ReMatch>,
    pub status: CloneStatus, // todo: Should this be an option?
    pub insert_loc: usize,
    /// Data for the insert. For example, used to draw its linear sequence.
    pub data_insert: Option<GenericData>,
    // note: This one may make more sense as a UI var.
    pub remove_stop_codons: bool,
    /// Work-in-progress cloning product sequence.
    pub product_seq: Seq,
    pub product_primers: Vec<Primer>,
}

impl Default for CloningState {
    fn default() -> Self {
        Self {
            insert_loc: 1,
            backbone_selected: Default::default(),
            backbone: Default::default(),
            res_common: Default::default(),
            re_matches_vec_common: Default::default(),
            re_matches_insert_common: Default::default(),
            status: Default::default(),
            data_insert: Default::default(),
            // backbone_data: Default::default(),
            remove_stop_codons: Default::default(),
            product_seq: Default::default(),
            product_primers: Vec::new(),
        }
    }
}

impl CloningState {
    /// Run this whenever something relevant in the cloning state changes. (e.g. selected a new backbone,
    /// a new insert location etc.) It synchronizes cloning product, validation, and RE matches.
    pub fn sync(
        &mut self,
        // mut seq_insert: &[Nucleotide],
        seq_insert: &mut Seq,
        backbone_lib: &[Backbone],
        re_lib: &[RestrictionEnzyme],
    ) {
        println!("Syncing cloning state...");
        // let backbone = self.get_backbone(backbone_lib);

        // todo: DRy with `get_backbone`, due to borrow errors. :(
        let backbone = match self.backbone_selected {
            BackboneSelected::Library(i) => {
                if i >= backbone_lib.len() {
                    eprintln!("Invalid index in backbone lib");
                    None
                } else {
                    Some(&backbone_lib[i])
                }
            }
            BackboneSelected::Opened => match self.backbone.as_ref() {
                Some(i) => Some(i),
                None => None,
            },
        };

        if let Some(backbone) = backbone {
            // todo: Put back. Getting a hang when this is there, eg from clicking auto pick insert loc.
            // (
            //     self.res_common,
            //     self.re_matches_vec_common,
            //     self.re_matches_insert_common,
            // ) = find_re_candidates(&backbone, seq_insert, &re_lib);

            self.product_seq = backbone.seq.clone();

            if seq_insert.len() < 3 {
                return;
            }

            // Trim off the stop codon at the end of the insert, A/R. This is quite specific.
            let insert_len = seq_insert.len();
            if self.remove_stop_codons
                && AminoAcid::from_codons(seq_insert[insert_len - 3..].try_into().unwrap())
                    == CodingResult::StopCodon
            {
                *seq_insert = seq_insert[..insert_len - 3].try_into().unwrap();
            }

            insert_into_seq(&mut self.product_seq, seq_insert, self.insert_loc).ok();
            self.status = CloneStatus::new(
                &backbone,
                self.insert_loc,
                seq_insert.len(),
                &self.product_seq,
            );
        }
    }

    pub fn get_backbone<'a>(&'a self, lib: &'a [Backbone]) -> Option<&'a Backbone> {
        match self.backbone_selected {
            BackboneSelected::Library(i) => {
                if i >= lib.len() {
                    eprintln!("Invalid index in backbone lib");
                    None
                } else {
                    Some(&lib[i])
                }
            }
            BackboneSelected::Opened => Some(self.backbone.as_ref().unwrap()),
        }
    }
}

#[derive(Clone, Copy, PartialEq)]
pub enum BackboneSelected {
    /// A backbone from our library; index.
    Library(usize),
    /// The current tab.
    Opened,
    // None, // todo: Rem this?
}

impl Default for BackboneSelected {
    fn default() -> Self {
        // Self::None
        Self::Opened
    }
}

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
    pub show_insert_picker: bool,
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
            state.ui.cloning_insert.seq_input = seq_to_str_lower(seq_this_ft);
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

    // state
    //     .ion_concentrations = state.ion_concentrations[state.active].clone());
    // state.file_active = None;
    state.portions.push(Default::default());
    state.volatile.push(Default::default());
    state.tabs_open.push(Default::default());
    state.ab1_data.push(Default::default());

    state.active = state.generic.len() - 1;

    // Make sure to create cloning primers before performing the insert, or the result will be wrong.
    make_cloning_primers(state);

    // todo: Unecessary clone? Due to borrow rules.
    let mut insert = state.ui.cloning_insert.seq_insert.clone();

    state.insert_nucleotides(&insert, state.cloning.insert_loc);

    let label = match state.ui.cloning_insert.feature_selected {
        Some(i) => state.ui.cloning_insert.features_loaded[i].label.clone(),
        None => "Cloning insert".to_owned(),
    };

    // todo: Eventually, we'll likely be pulling in sequences already associated with a feature;
    // todo: Use the already existing data instead.
    state.generic[state.active].features.push(Feature {
        range: RangeIncl::new(
            state.cloning.insert_loc,
            state.cloning.insert_loc + insert.len() - 1,
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

/// Check if the (eg His) tag is in frame with the start of the coding region, and that there is
/// no stop codon between the end of the coding region, and start of the sequence.
///
/// We assume, for now, that the insert location is the start of a coding region.
/// Note: `his_tag` here must have already been shifted based on the insert.
/// todo: This assumption may not be always valid!
// fn tag_in_frame(backbone: &Backbone, coding_region_start: usize, insert_len: usize) -> Status {
fn tag_in_frame(
    seq_product: &[Nucleotide],
    his_tag: &Option<RangeIncl>,
    coding_region_start: usize,
) -> Status {
    if seq_product.is_empty() {
        eprintln!("Error checking tag frame: Product sequence is empty");
        return Status::NotApplicable;
    }
    match his_tag {
        Some(tag_range) => {
            let tag_start = tag_range.start;
            let cr_start = coding_region_start % seq_product.len();

            //
            // let (tag_start, cr_start) = match backbone.direction {
            //     PrimerDirection::Forward => {
            //
            //         (tag_range.start + insert_len, cr_start)
            //     }
            //     // todo: QC this.
            //     PrimerDirection::Reverse => {
            //         let tag_end = tag_range.end % backbone.seq.len();
            //         (coding_region_start, tag_end + insert_len)
            //     }
            // };

            if cr_start > tag_start {
                eprintln!("Error with insert loc and tag end. Coding region start: {cr_start}, tag start: {tag_start}");
                return Status::Fail;
            }

            println!("Tag: {tag_start}, CR: {coding_region_start}");

            // todo: Helper fn for this sort of check in general?
            // If there is a stop codon between the end of the sequence and the
            for i in (cr_start..tag_start).step_by(3) {
                if i + 2 >= seq_product.len() {
                    eprintln!("Error: Invalid backbone seq in his check");
                    return Status::Fail;
                }

                if let CodingResult::StopCodon = AminoAcid::from_codons([
                    // Offset for 1-based indexing.
                    seq_product[i - 1],
                    seq_product[i + 0],
                    seq_product[i + 1],
                ]) {
                    // todo: Confirm this is the full seq with insert adn vec.
                    println!("Stop codon between seq start and his tag found at index {i}"); // todo temp
                    return Status::Fail;
                }
            }

            let dist_from_start_codon = tag_start - cr_start;

            if dist_from_start_codon % 3 == 0 {
                Status::Pass
            } else {
                Status::Fail
            }
        }
        None => Status::NotApplicable,
    }
}

/// For a given insert and vector, find suitable  restriction enzymes for cloning.
/// Make sure that the insert sequence is properly buffered to allow for RE matching outside
/// of the coding region (etc)'s range, upstream of this.
/// todo: Fow now, we limit our options to unique-match, single-cutters.
///
/// Returns REs, matches vector, matches insert.
pub fn find_re_candidates<'a>(
    backbone: &Backbone,
    seq_insert: &[Nucleotide],
    // re_match_sets: &[&[ReMatch]],
    // insert_tab: usize,
    // insert_range: RangeIncl,
    lib: &'a [RestrictionEnzyme],
    // ) -> Vec<&'a RestrictionEnzyme> {
) -> (Vec<RestrictionEnzyme>, Vec<ReMatch>, Vec<ReMatch>) {
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

    let res = result.into_iter().map(|r| r.clone()).collect();

    let mut matches_vector_common = Vec::new();
    let mut matches_insert_common = Vec::new();

    for match_insert in &matches_insert {
        for match_bb in &matches_backbone {
            if match_insert.lib_index == match_bb.lib_index {
                matches_vector_common.push(match_bb.clone());
                matches_insert_common.push(match_insert.clone());
            }
        }
    }

    (res, matches_vector_common, matches_insert_common)
}

/// Validation checks for cloning.
#[derive(Default)]
pub struct CloneStatus {
    pub rbs_dist: Status,
    pub downstream_of_promoter: Status,
    pub upstream_of_terminator: Status,
    // todo: Check for terminator?
    pub direction: Status,
    pub tag_frame: Status,
    // pub primer_quality: Status,
    // pub re_dist: Status,
}

impl CloneStatus {
    pub fn new(
        backbone: &Backbone,
        insert_loc: usize,
        insert_len: usize,
        seq_product: &[Nucleotide],
    ) -> Self {
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
                if insert_loc < p.end {
                    Status::Pass
                } else {
                    Status::Fail
                }
            }
            None => Status::NotApplicable,
        };

        let direction = Status::Pass; // todo

        // let tag_frame = tag_in_frame(&backbone, insert_loc, insert_len);
        // Note: We are assuming for now, that the coding region start/frame, is the insert location.

        let his_tag_shifted = match backbone.his_tag {
            Some(tag) => Some(RangeIncl::new(tag.start + insert_len, tag.end + insert_len)),
            None => None,
        };

        let tag_frame = tag_in_frame(seq_product, &his_tag_shifted, insert_loc);

        Self {
            rbs_dist,
            downstream_of_promoter,
            upstream_of_terminator,
            direction,
            tag_frame,
        }
    }
}
