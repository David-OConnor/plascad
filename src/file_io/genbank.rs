//! Parse and write GenBank files. This converts between this format for sequences, features,
//! primers, notes etc., and our own equivalents. We use the [gb_io](https://docs.rs/gb-io/latest/gb_io)
//! library.
//!
//! (NIH article on GenBank)[https://www.ncbi.nlm.nih.gov/genbank/]

use std::{
    collections::HashMap,
    fs::{File, OpenOptions},
    io::{self, ErrorKind, Write},
    ops::Range,
    path::Path,
};

use gb_io::{
    self,
    reader::SeqReader,
    seq::{After, Before, Location},
    writer::SeqWriter,
};

use crate::{
    file_io::GenericData,
    primer::{Primer, PrimerDirection},
    sequence::{
        seq_complement, Feature, FeatureDirection, FeatureType, Nucleotide,
        SeqTopology,
    },
    Metadata, Reference,
};
use crate::primer::PrimerData;

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
        let mut primers = Vec::new();

        // This is almost awkward enough to write a parser instead of using gb_io.
        for feature in &seq.features {
            let feature_type = FeatureType::from_external_str(&feature.kind.to_string());

            // We parse label from qualifiers.
            // I'm unsure how direction works in GenBank files. It appears it's some mix of the LEFT/RIGHT
            // qualifiers, feature type, and if the location is complement, or forward.
            let mut direction = FeatureDirection::None;
            let mut label = String::new();

            let index_range = match &feature.location {
                // gb_io seems to list the start of the range as 1 too early; compensate.
                Location::Range(start, end) => (start.0 as usize + 1, end.0 as usize),
                Location::Complement(inner) => match **inner {
                    Location::Range(start, end) => {
                        direction = FeatureDirection::Reverse;
                        (start.0 as usize + 1, end.0 as usize)
                    }
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

            for v in feature.qualifier_values("label".into()) {
                label = v.to_owned();
                break;
            }

            // Todo: Update or remove this A/R.
            match feature_type {
                FeatureType::Primer => {
                    if direction != FeatureDirection::Reverse {
                        direction = FeatureDirection::Forward;
                    }
                }
                FeatureType::CodingRegion => direction = FeatureDirection::Forward,
                _ => (),
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

            // GenBank stores primer bind sites (Which we treat as volatile), vice primer sequences.
            // Infer the sequence using the bind indices, and the main sequence.
            if feature_type == FeatureType::Primer {
                let sequence = match direction {
                    FeatureDirection::Reverse => {
                        let compl = seq_complement(&seq_);
                        compl[seq_.len() - (index_range.1 - 1)..seq_.len() - (index_range.0)]
                            .to_vec()
                    }
                    // See other notes on start range index being odd.
                    _ => seq_[index_range.0 - 1..index_range.1].to_vec(),
                };

                let volatile = PrimerData::new(&sequence);
                primers.push(Primer {
                    sequence,
                    name: label,
                    description: None, // todo: Can we populate this?
                    volatile,
                });
                continue;
            }

            // Parse notes from qualifiers other than label and direction.
            let mut notes = HashMap::new();
            for (qual_key, val) in &feature.qualifiers {
                if qual_key == "label" {
                    continue; // We handle this separately.
                }
                if let Some(v) = val {
                    notes.insert(qual_key.to_string(), v.clone());
                }
            }

            features.push(Feature {
                index_range,
                feature_type,
                direction,
                label,
                color_override: None,
                notes,
            })
        }

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

        let (source, organism) = match seq.source {
            Some(src) => (Some(src.source), src.organism),
            None => (None, None),
        };

        let metadata = Metadata {
            plasmid_name: seq.keywords.clone().unwrap_or(String::new()),
            comments: seq.comments.clone(),
            definition: seq.definition.clone(),
            accession: seq.accession.clone(),
            version: seq.version.clone(),
            keywords: seq.keywords.clone(),
            locus: seq.name.clone().unwrap_or_default(),
            source,
            organism,
            references,
            ..Default::default()
        };

        return Ok(GenericData {
            seq: seq_,
            topology,
            features,
            primers,
            metadata,
        });
    }

    Err(io::Error::new(
        ErrorKind::InvalidData,
        format!("No GenBank sequences found"),
    ))
}

/// Export our local state into the GenBank format. This includes sequence, features, and primers.
pub fn export_genbank(
    data: &GenericData,
    primer_matches: &[(PrimerDirection, Range<usize>, String)],
    path: &Path,
) -> io::Result<()> {
    let file = OpenOptions::new()
        .write(true)
        .create(true) // Create the file if it doesn't exist
        .open(path)?;

    let mut gb_data = gb_io::seq::Seq::empty();

    gb_data.seq = data.seq.iter().map(|nt| nt.to_u8_letter()).collect();

    gb_data.topology = match data.topology {
        SeqTopology::Circular => gb_io::seq::Topology::Circular,
        SeqTopology::Linear => gb_io::seq::Topology::Linear,
    };

    for feature in &data.features {
        let mut qualifiers = vec![("label".into(), Some(feature.label.clone()))];

        for note in &feature.notes {
            qualifiers.push(((&**note.0).into(), Some(note.1.clone())));
        }

        match feature.direction {
            FeatureDirection::Forward => {
                qualifiers.push(("direction".into(), Some("right".to_owned())))
            }
            FeatureDirection::Reverse => {
                qualifiers.push(("direction".into(), Some("left".to_owned())))
            }
            _ => (),
        }

        let start: i64 = feature.index_range.0.try_into().unwrap();
        let end: i64 = feature.index_range.1.try_into().unwrap();

        let location = match feature.direction {
            FeatureDirection::Reverse => Location::Complement(Box::new(Location::Range(
                (end + 1, Before(false)),
                (start, After(false)), // todo: Offset on end. Whhy?
            ))),
            _ => Location::Range(
                (start - 1, Before(false)),
                (feature.index_range.1.try_into().unwrap(), After(false)),
            ),
        };

        gb_data.features.push(gb_io::seq::Feature {
            kind: feature.feature_type.to_external_str().into(),
            location,
            qualifiers,
        });
    }

    for (dir, indexes, name) in primer_matches {
        // todo: Location code is DRY with features.
        let start: i64 = indexes.start.try_into().unwrap();
        let end: i64 = indexes.end.try_into().unwrap();

        let location = match dir {
            PrimerDirection::Forward => Location::Range(
                (start, Before(false)),
                (indexes.end.try_into().unwrap(), After(false)),
            ),
            PrimerDirection::Reverse => Location::Complement(Box::new(Location::Range(
                (data.seq.len() as i64 - (end + 1), Before(false)),
                (data.seq.len() as i64 - (start - 1), After(false)), // todo: Offset on end. Whhy?
            ))),
        };

        gb_data.features.push(gb_io::seq::Feature {
            kind: "primer_bind".into(),
            location,
            qualifiers: vec![("label".into(), name.to_owned().into())],
        });
    }

    let md = &data.metadata;

    gb_data.comments = md.comments.clone();
    gb_data.source = Some(gb_io::seq::Source {
        source: md.source.clone().unwrap_or_default(),
        organism: md.organism.clone(),
    });

    if md.source.is_none() && md.organism.is_none() {
        gb_data.source = None;
    }

    // data.keywords = md.keywords.clone();
    gb_data.keywords = Some(md.plasmid_name.clone());
    gb_data.version = md.version.clone();
    gb_data.accession = md.accession.clone();
    gb_data.definition = md.definition.clone();
    gb_data.name = Some(md.locus.clone());

    for ref_ in &md.references {
        gb_data.references.push(gb_io::seq::Reference {
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
    writer.write(&gb_data)
}
