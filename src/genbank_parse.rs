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

use gb_io::{self, reader::SeqReader};

use crate::{
    primer::{Primer, PrimerData},
    sequence::{seq_to_str, Feature, FeatureDirection, FeatureType, Nucleotide, Seq, SeqTopology},
};

#[derive(Default)]
pub struct GenBankData {
    pub seq: Seq,
    pub topology: SeqTopology,
    pub features: Vec<Feature>,
    // pub primers: Option<Vec<Primer>>,
}

// fn export_features(buf: &mut Vec<u8>, features: &[Feature]) -> io::Result<()> {
//     Ok(())
// }
//
// fn export_primers(buf: &mut Vec<u8>, primers: &[PrimerData]) -> io::Result<()> {
//     Ok(())
// }

/// Read a file in the GenBank format.
/// [Rust docs ref of fields](https://docs.rs/gb-io/latest/gb_io/seq/struct.Seq.html)
pub fn import_genbank(path: &Path) -> io::Result<GenBankData> {
    let mut result = GenBankData::default();

    let mut file = File::open(path)?;

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

        // todo: No primers?
        // todo: What is contig location? Do we need that?

        return Ok(GenBankData {
            seq: seq_,
            topology,
            features,
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
    topology: SeqTopology,
    features: &[Feature],
    primers: &[PrimerData],
    path: &Path,
) -> io::Result<()> {
    let mut file = OpenOptions::new()
        .write(true)
        .create(true) // Create the file if it doesn't exist
        .open(path)?;

    let mut buf = Vec::new();

    buf.extend(seq_to_str(seq).as_bytes());

    // export_features(&mut buf, features)?;
    // export_primers(&mut buf, primers)?;

    file.write_all(&buf)?;

    Ok(())
}
