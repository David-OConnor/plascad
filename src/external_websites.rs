//! For opening the browser to NCBI BLAST, PDB etc.
//!
//! PDB Search API: https://search.rcsb.org/#search-api
//! PDB Data API: https://data.rcsb.org/#data-api

use std::{char::MAX, time::Duration};

use bincode::{Decode, Encode};
use reqwest::{self, header::CONTENT_TYPE};
use serde::{Deserialize, Serialize};
use serde_json::{self, Value};
use url::Url;

use crate::{
    protein::{aa_seq_to_str, Protein},
    sequence::{seq_to_str, Nucleotide},
    Selection, State,
};

const NCBI_BLAST_URL: &str = "https://blast.ncbi.nlm.nih.gov/Blast.cgi";

const PDB_BASE_URL: &str = "https://www.rcsb.org/structure";
const PDB_3D_VIEW_URL: &str = "https://www.rcsb.org/3d-view";
const PDB_STRUCTURE_FILE_URL: &str = "https://files.rcsb.org/view";

const PDB_SEARCH_API_URL: &str = "https://search.rcsb.org/rcsbsearch/v2/query";
const PDB_DATA_API_URL: &str = "https://data.rcsb.org/rest/v1/core/entry";

// An arbitrary limit to prevent excessive queries to the PDB data api,
// and to simplify display code.
const MAX_PDB_RESULTS: usize = 8;

/// BLAST the selected Feature, primer, or selection. Prioritize the selection.
/// This function handles extracting the sequence to BLAST from possible selections.
pub fn blast(state: &State) {
    let val = match state.ui.text_selection {
        Some(sel) => {
            // Don't format sel directly, as we insert the bp count downstream for use with feature selections.
            Some((
                sel.index_seq(&state.generic.seq),
                format!(
                    "{}, {}..{}",
                    state.generic.metadata.plasmid_name, sel.start, sel.end
                ),
            ))
        }
        None => match state.ui.selected_item {
            Selection::Feature(feat_i) => {
                if state.generic.features.len() < feat_i + 1 {
                    eprintln!("Invalid selected feature");
                    None
                } else {
                    let feature = &state.generic.features[feat_i];
                    Some((
                        feature.range.index_seq(&state.generic.seq),
                        feature.label.clone(),
                    ))
                }
            }
            Selection::Primer(prim_i) => {
                if state.generic.primers.len() < prim_i + 1 {
                    eprintln!("Invalid selected primer");
                    None
                } else {
                    let primer = &state.generic.primers[prim_i];
                    Some((Some(&primer.sequence[..]), primer.name.clone()))
                }
            }
            Selection::None => None,
        },
    };

    // todo: Handle reverse.

    if let Some((seq, name)) = val {
        if let Some(s) = seq {
            open_blast(s, &name);
        }
    }
}

/// Open the web browser to a NCBI-BLAST page, of the sequence of interest.
///
///Example BLAST
/// note: There appears to be a NT limit that will fail most full plastmids when using the GET api.
/// ?PAGE_TYPE=BlastSearch&CMD=Web&LAYOUT=OneWindow&PROGRAM=blastn&MEGABLAST=on&PAGE=Nucleotides&DATABASE=nr
/// &FORMAT_TYPE=HTML&NCBI_GI=on&SHOW_OVERVIEW=on&QUERY=%3Ettt%20%20(43%20..%20905%20%3D%20863%20bp)
/// %0AACTCACTATAGGGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACCGGTATGACTAGTATGGAAGACGCCAAAAACATAAAGAAAGGCCCG
/// GCGCCATTCTATCCGCTGGAAGATGGAACCGCTGGAGAGCAACTGCATAAGGCTATGAAGAGATACGCCCTGGTTCCTGGAACAATTGCTTTTACAGA
/// TGCACATATCGAGGTGGACATCACTTACGCTGAGTACTTCGAAATGTCCGTTCGGTTGGCAGAAGCTATGAAACGATATGGGCTGAATACAAATCACAGA
/// ATCGTCGTATGCAGTGAAAACTCTCTTCAATTCTTTATGCCGGTGTTGGGCGCGTTATTTATCGGAGTTGCAGTTGCGCCCGCGAACGACATTTATAATGA
/// ACGTGAATTGCTCAACAGTATGGGCATTTCGCAGCCTACCGTGGTGTTCGTTTCCAAAAAGGGGTTGCAAAAAATTTTGAACGTGCAAAAAAAGCTCCCAAT
/// CATCCAAAAAATTATTATCATGGATTCTAAAACGGATTACCAGGGATTTCAGTCGATGTACACGTTCGTCACATCTCATCTACCTCCCGGTTTTAATGAATAC
/// GATTTTGTGCCAGAGTCCTTCGATAGGGACAAGACAATTGCACTGATCATGAACTCCTCTGGATCTACTGGTCTGCCTAAAGGTGTCGCTCTGCCTCATAGAACT
/// GCCTGCGTGAGATTCTCGCATGCCAGAGATCCTATTTTTGGCAATCAAATCATTCCGGATACTGCGATTTTAAGTGTTGTTCCATTCCATCACGGTTTTGGAA
/// TGTTTACTACACTCGGATATTTGATATGTGGATTTCGAGTCGTCTTAATGTATAGAT
fn open_blast(seq: &[Nucleotide], seq_name: &str) {
    let text_query = format!(">{seq_name} ({} bp)\n{}", seq.len(), seq_to_str(seq));

    let params = vec![
        ("PAGE_TYPE", "BlastSearch"),
        ("CMD", "Web"),
        ("LAYOUT", "OneWindow"),
        ("PROGRAM", "blastn"),
        ("MEGABLAST", "on"),
        ("PAGE", "Nucleotides"),
        ("DATABASE", "nr"),
        ("FORMAT_TYPE", "HTML"),
        ("NCBI_GI", "on"),
        ("SHOW_OVERVIEW", "on"),
        ("QUERY", &text_query),
    ];

    let mut url = Url::parse(NCBI_BLAST_URL).unwrap();
    {
        let mut query_pairs = url.query_pairs_mut();
        for (key, value) in params {
            query_pairs.append_pair(key, value);
        }
    }

    // Open the URL in the default web browser
    if let Err(e) = webbrowser::open(url.as_str()) {
        eprintln!("Failed to open the web browser: {:?}", e);
    }
}

#[derive(Default, Serialize)]
struct PdbSearchParams {
    value: String,
    sequence_type: String,
    evalue_cutoff: u8,
    identity_cutoff: f32,
}

#[derive(Default, Serialize)]
struct PdbSearchQuery {
    #[serde(rename = "type")]
    type_: String,
    service: String,
    parameters: PdbSearchParams,
}

#[derive(Default, Serialize)]
struct SearchRequestOptions {
    scoring_strategy: String,
}

#[derive(Default, Serialize)]
struct PdbPayloadSearch {
    return_type: String,
    query: PdbSearchQuery,
    #[serde(skip_serializing_if = "Option::is_none")]
    request_options: Option<SearchRequestOptions>,
    #[serde(skip_serializing_if = "Option::is_none")]
    request_info: Option<String>,
}

#[derive(Default, Debug, Encode, Decode, Deserialize)]
pub struct PdbSearchResult {
    pub identifier: String,
    pub score: f32,
}

#[derive(Default, Debug, Deserialize)]
struct PdbSearchResults {
    query_id: String,
    result_type: String,
    total_count: u32,
    result_set: Vec<PdbSearchResult>,
}

#[derive(Default, Debug, Deserialize)]
struct PdbStruct {
    title: String,
}

#[derive(Default, Debug, Deserialize)]
struct PdbDataResults {
    #[serde(rename = "struct")]
    struct_: PdbStruct,
}

/// This doesn't deserialize directly; it's the format we use internally.
#[derive(Encode, Decode)]
pub struct PdbData {
    pub rcsb_id: String,
    pub title: String,
}

/// Load PDB data using [its API](https://search.rcsRb.org/#search-api)
/// Returns the set of PDB ID matches, with scores.
pub fn load_pdb_data(protein: &Protein) -> Result<Vec<PdbData>, reqwest::Error> {
    let payload_search = PdbPayloadSearch {
        return_type: "entry".to_string(),
        query: PdbSearchQuery {
            type_: "terminal".to_owned(),
            service: "sequence".to_owned(),
            parameters: PdbSearchParams {
                value: aa_seq_to_str(&protein.aa_seq),
                sequence_type: "protein".to_owned(),
                evalue_cutoff: 1,
                identity_cutoff: 0.9,
            },
        },
        request_options: Some(SearchRequestOptions {
            scoring_strategy: "sequence".to_owned(),
        }),

        // return_type: "assembly".to_string(), // todo: Experiment.
        ..Default::default()
    };

    // todo: Limit the query to our result cap, instead of indexing after?

    let payload_json = serde_json::to_string(&payload_search).unwrap();

    let client = reqwest::blocking::Client::builder()
        .timeout(Duration::from_secs(4))
        .build()?;

    // todo: New thread for this A/R.

    let resp = client
        .post(PDB_SEARCH_API_URL)
        .header(CONTENT_TYPE, "application/json")
        .body(payload_json)
        .send()?;

    let search_data: PdbSearchResults = resp.json()?;

    let mut result_search = Vec::new();
    for (i, r) in search_data.result_set.into_iter().enumerate() {
        if i < MAX_PDB_RESULTS {
            result_search.push(r);
        }
    }

    let mut result = Vec::new();
    for r in result_search {
        let resp = client
            .get(&format!("{PDB_DATA_API_URL}/{}", r.identifier))
            // .header(CONTENT_TYPE, "application/json")
            .send()?;

        // println!("resp data: {:?}", resp.text());
        let data: PdbDataResults = resp.json()?;
        // println!("PDB Data: {:?}", data);

        result.push(PdbData {
            rcsb_id: r.identifier,
            title: data.struct_.title,
        })
    }

    // Now, load data for each result.t

    Ok(result)
}

/// Open a PDB search for this protein's sequence, given a PDB ID, which we load from the API.
pub fn open_pdb(pdb_id: &str) {
    // Open the URL in the default web browser
    if let Err(e) = webbrowser::open(&format!("{PDB_BASE_URL}/{pdb_id}")) {
        eprintln!("Failed to open the web browser: {:?}", e);
    }
}

/// Open a PDB search for this protein's sequence, given a PDB ID, which we load from the API.
pub fn open_pdb_3d_view(pdb_id: &str) {
    // Open the URL in the default web browser
    if let Err(e) = webbrowser::open(&format!("{PDB_3D_VIEW_URL}/{pdb_id}")) {
        eprintln!("Failed to open the web browser: {:?}", e);
    }
}

/// Load PDB structure data in the PDBx/mmCIF format. This is a modern, text-based format.
/// It avoids the XML, and limitations of the other two available formats.
/// todo: When to use wwpdb vs rscb?
pub fn load_pdb_structure(pdb_id: &str) {
    let pdb_id = pdb_id.to_owned().to_lowercase();

    // todo: Use EGUI_file to save the file and reqwest to load it.

    // let url_pdb_format = format!("https://files.wwpdb.org/pub/pdb/data/structures/divided/pdb/zg/pdb{pdb_id}.ent.gz");
    // let url_pdbx_format = format!("https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/zg/{pdb_id}.cif.gz");

    let url = format!("{PDB_STRUCTURE_FILE_URL}/{pdb_id}.cif");

    if let Err(e) = webbrowser::open(&url) {
        eprintln!("Failed to open the web browser: {:?}", e);
    }
}
