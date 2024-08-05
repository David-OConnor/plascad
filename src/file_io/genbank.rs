//! Parse and write GenBank files. This converts between this format for sequences, features,
//! primers, notes etc., and our own equivalents. We use the [gb_io](https://docs.rs/gb-io/latest/gb_io)
//! library.
//!
//! (NIH article on GenBank)[https://www.ncbi.nlm.nih.gov/genbank/]

use std::{
    collections::HashMap,
    fs::{File, OpenOptions},
    io::{self, ErrorKind, Write},
    path::Path,
};

use gb_io::{self, reader::SeqReader, writer::SeqWriter};

use crate::{
    file_io::GenericData,
    primer::{Primer, PrimerData},
    sequence::{seq_to_str, Feature, FeatureDirection, FeatureType, Nucleotide, Seq, SeqTopology},
    Reference,
};
// fn export_features(buf: &mut Vec<u8>, features: &[Feature]) -> file_io::Result<()> {
//     Ok(())
// }
//
// fn export_primers(buf: &mut Vec<u8>, primers: &[PrimerData]) -> file_io::Result<()> {
//     Ok(())
// }

/// Read a file in the GenBank format.
/// [Rust docs ref of fields](https://docs.rs/gb-io/latest/gb_io/seq/struct.Seq.html)
pub fn import_genbank(path: &Path) -> io::Result<GenericData> {
    let file = File::open(path)?;

    // todo: This currently only handles a single sequene. It returns the first found.

    // todo: Should we take advantage of other features?
    //     pub name: Option<String>,
    //     pub topology: Topology,
    //     pub date: Option<Date>,
    //     pub len: Option<usize>,
    //     pub molecule_type: Option<String>,
    //     pub division: String,
    //     pub definition: Option<String>,
    //     pub accession: Option<String>,
    //     pub version: Option<String>,
    //     pub source: Option<Source>,
    //     pub dblink: Option<String>,
    //     pub keywords: Option<String>,
    //     pub references: Vec<Reference>,
    //     pub comments: Vec<String>,
    //     pub seq: Vec<u8>,
    //     pub contig: Option<Location>,
    //     pub features: Vec<Feature>,
    // }

    for seq in SeqReader::new(file) {
        let seq = seq.map_err(|e| {
            io::Error::new(
                ErrorKind::InvalidData,
                format!("Unable to get GenBank seq {e}"),
            )
        })?;

        let mut seq_ = Vec::new();

        for nt in &seq.seq {
            match Nucleotide::from_u8_letter(*nt) {
                Ok(n) => seq_.push(n),
                Err(_) => {
                    eprintln!("Unexpected char in GenBank sequence: {:?}", nt);
                }
            }
        }

        let topology = match seq.topology {
            gb_io::seq::Topology::Linear => SeqTopology::Linear,
            gb_io::seq::Topology::Circular => SeqTopology::Circular,
        };

        let mut features = Vec::new();

        // This is almost awkward enough to write a parser instead of using gb_io.
        for feature in &seq.features {
            // We parse label and direction from qualifiers.
            let mut direction = FeatureDirection::None;
            let mut label = String::new();

            let index_range = match &feature.location {
                gb_io::seq::Location::Range(start, end) => (start.0 as usize, end.0 as usize),
                gb_io::seq::Location::Complement(inner) => match **inner {
                    gb_io::seq::Location::Range(start, end) => (start.0 as usize, end.0 as usize),
                    _ => {
                        eprintln!("Unexpected gb_io compl range type: {:?}", feature.location);
                        (0, 0)
                    }
                },
                _ => {
                    eprintln!("Unexpected gb_io range type: {:?}", feature.location);
                    (0, 0)
                }
            };

            // These Atom and sets from gb_io are tough to work with.
            for v in feature.qualifier_values("label".into()) {
                label = v.to_owned();
                break;
            }

            for v in feature.qualifier_values("direction".into()) {
                let v = v.to_lowercase();
                if v == "right" {
                    direction = FeatureDirection::Forward;
                    break;
                } else if v == "left" {
                    direction = FeatureDirection::Reverse;
                    break;
                }
            }

            // Parse notes from qualifiers other than label and direction.
            let mut notes = HashMap::new();
            for qual in &feature.qualifiers {}

            // todo: Handle primers from the `primer_bind` feature kind.

            features.push(Feature {
                index_range,
                feature_type: FeatureType::from_external_str(&feature.kind.to_string()),
                direction,
                label,
                color_override: None,
                notes,
            })
        }

        let mut plasmid_name = String::new(); // todo
        let mut primers = Vec::new(); // todo

        let mut references = Vec::new();
        for ref_ in &seq.references {
            references.push(Reference {
                description: ref_.description.clone(),
                authors: ref_.authors.clone(),
                consortium: ref_.consortium.clone(),
                title: ref_.title.clone(),
                journal: ref_.journal.clone(),
                pubmed: ref_.pubmed.clone(),
                remark: ref_.remark.clone(),
            })
        }

        // todo: No primers?
        // todo: What is contig location? Do we need that?

        return Ok(GenericData {
            seq: seq_,
            plasmid_name,
            topology,
            features,
            primers,
            comments: seq.comments.clone(),
            references,
        });
    }

    Err(io::Error::new(
        ErrorKind::InvalidData,
        format!("No GenBank sequences found"),
    ))
}

/// Export our local state into the GenBank format. This includes sequence, features, and primers.
pub fn export_genbank(
    seq: &[Nucleotide],
    plasmid_name: &str,
    topology: SeqTopology,
    features: &[Feature],
    primers: &[PrimerData],
    comments: &[String],
    references: &[Reference],
    path: &Path,
) -> io::Result<()> {
    let file = OpenOptions::new()
        .write(true)
        .create(true) // Create the file if it doesn't exist
        .open(path)?;

    let mut data = gb_io::seq::Seq::empty();

    data.seq = seq.iter().map(|nt| nt.to_u8_letter()).collect();

    data.topology = match topology {
        SeqTopology::Circular => gb_io::seq::Topology::Circular,
        SeqTopology::Linear => gb_io::seq::Topology::Linear,
    };

    // todo: Features
    // todo: Handle primers.

    data.comments = comments.to_vec().clone();

    for ref_ in references {
        data.references.push(gb_io::seq::Reference {
            description: ref_.description.clone(),
            authors: ref_.authors.clone(),
            consortium: ref_.consortium.clone(),
            title: ref_.title.clone(),
            journal: ref_.journal.clone(),
            pubmed: ref_.pubmed.clone(),
            remark: ref_.remark.clone(),
        })
    }

    let mut writer = SeqWriter::new(file);
    writer.write(&data)
}
