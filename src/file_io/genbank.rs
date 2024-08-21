//! Parse and write GenBank files. This converts between this format for sequences, features,
//! primers, notes etc., and our own equivalents. We use the [gb_io](https://docs.rs/gb-io/latest/gb_io)
//! library.
//!
//! (NIH article on GenBank)[https://www.ncbi.nlm.nih.gov/genbank/]

use std::{
    fs::File,
    io::{self, ErrorKind},
    path::Path,
};

use gb_io::{
    self,
    reader::SeqReader,
    seq::{After, Before, Location},
    writer::SeqWriter,
};

use crate::{
    file_io::{get_filename, GenericData},
    primer::{Primer, PrimerData, PrimerDirection, PrimerMatch},
    sequence::{
        seq_complement, Feature, FeatureDirection, FeatureType, Metadata, Nucleotide, Reference,
        SeqTopology,
    },
    util::RangeIncl,
};

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

        let (features, primers) = parse_features_primers(&seq.features, &seq_);

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

        // todo: What is contig location? Do we need that?

        let (source, organism) = match seq.source {
            Some(src) => (Some(src.source), src.organism),
            None => (None, None),
        };

        let metadata = Metadata {
            plasmid_name: get_filename(path),
            comments: seq.comments.clone(),
            definition: seq.definition.clone(),
            accession: seq.accession.clone(),
            version: seq.version.clone(),
            keywords: seq.keywords.clone(),
            locus: seq.name.clone().unwrap_or_default(),
            source,
            organism,
            references,
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
        "No GenBank sequences found",
    ))
}

/// Parse features and primers, from GenBank's feature list.
fn parse_features_primers(
    features: &[gb_io::seq::Feature],
    seq: &[Nucleotide],
) -> (Vec<Feature>, Vec<Primer>) {
    let mut result_ft = Vec::new();
    let mut primers = Vec::new();

    let compl = seq_complement(seq);

    for feature in features {
        let feature_type = FeatureType::from_external_str(feature.kind.as_ref());

        // We parse label from qualifiers.
        // I'm unsure how direction works in GenBank files. It appears it's some mix of the LEFT/RIGHT
        // qualifiers, feature type, and if the location is complement, or forward.
        let mut direction = FeatureDirection::None;
        let mut label = String::new();

        let mut range = match &feature.location {
            // gb_io seems to list the start of the range as 1 too early; compensate.
            Location::Range(start, end) => RangeIncl::new(start.0 as usize + 1, end.0 as usize),
            Location::Complement(inner) => match **inner {
                Location::Range(start, end) => {
                    direction = FeatureDirection::Reverse;
                    RangeIncl::new(start.0 as usize + 1, end.0 as usize)
                }
                _ => {
                    eprintln!("Unexpected gb_io compl range type: {:?}", feature.location);
                    RangeIncl::new(1, 1)
                }
            },
            _ => {
                eprintln!("Unexpected gb_io range type: {:?}", feature.location);
                RangeIncl::new(1, 1)
            }
        };

        for v in feature.qualifier_values("label".into()) {
            v.clone_into(&mut label);
            break;
        }

        match feature_type {
            FeatureType::Primer => {
                if direction != FeatureDirection::Reverse {
                    direction = FeatureDirection::Forward;
                }
            }
            FeatureType::CodingRegion => {
                // As CDS regions are always directional, if not identified as reverse due to the range
                // being in complement format, set it to forward.
                if direction == FeatureDirection::None {
                    direction = FeatureDirection::Forward
                }
            }
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

        // Passing through the origin, most likely. Loop around past the end.
        if range.end < range.start {
            range = RangeIncl::new(range.start, seq.len() + range.end);
        }

        // GenBank stores primer bind sites (Which we treat as volatile), vice primer sequences.
        // Infer the sequence using the bind indices, and the main sequence.
        if feature_type == FeatureType::Primer {
            let sequence = match direction {
                FeatureDirection::Reverse => {
                    let range =
                        RangeIncl::new(seq.len() - (range.end) + 1, seq.len() - (range.start));

                    range.index_seq(&compl).unwrap_or_default()
                }
                _ => range.index_seq(&seq).unwrap_or_default(),
            }
            .to_vec();

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
        // let mut notes = HashMap::new();
        let mut notes = Vec::new();
        for (qual_key, val) in &feature.qualifiers {
            if qual_key == "label" {
                continue; // We handle this separately.
            }
            if let Some(v) = val {
                // notes.insert(qual_key.to_string(), v.clone());
                notes.push((qual_key.to_string(), v.clone()));
            }
        }

        result_ft.push(Feature {
            range,
            feature_type,
            direction,
            label,
            color_override: None,
            notes,
        })
    }

    (result_ft, primers)
}

/// Export our local state into the GenBank format. This includes sequence, features, and primers.
pub fn export_genbank(
    data: &GenericData,
    primer_matches: &[(PrimerMatch, String)],
    path: &Path,
) -> io::Result<()> {
    let file = File::create(path)?;

    let mut gb_data = gb_io::seq::Seq::empty();

    gb_data.seq = data.seq.iter().map(|nt| nt.to_u8_letter()).collect();

    gb_data.topology = match data.topology {
        SeqTopology::Circular => gb_io::seq::Topology::Circular,
        SeqTopology::Linear => gb_io::seq::Topology::Linear,
    };

    for feature in &data.features {
        let mut qualifiers = vec![("label".into(), Some(feature.label.clone()))];

        for note in &feature.notes {
            qualifiers.push(((&*note.0).into(), Some(note.1.clone())));
            // qualifiers.push((note.0.clone(), Some(note.1.clone())));
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

        let start: i64 = feature.range.start.try_into().unwrap();
        let end: i64 = feature.range.end.try_into().unwrap();

        let location = match feature.direction {
            FeatureDirection::Reverse => Location::Complement(Box::new(Location::Range(
                (start - 1, Before(false)),
                (end, After(false)),
            ))),
            _ => Location::Range(
                (start - 1, Before(false)),
                (feature.range.end.try_into().unwrap(), After(false)),
            ),
        };

        gb_data.features.push(gb_io::seq::Feature {
            kind: feature.feature_type.to_external_str().into(),
            location,
            qualifiers,
        });
    }

    for (prim_match, name) in primer_matches {
        // todo: Location code is DRY with features.
        let start: i64 = prim_match.range.start.try_into().unwrap();

        let location = match prim_match.direction {
            PrimerDirection::Forward => Location::Range(
                (start - 1, Before(false)),
                (prim_match.range.end.try_into().unwrap(), After(false)),
            ),
            PrimerDirection::Reverse => Location::Complement(Box::new(Location::Range(
                (start - 1, Before(false)),
                (prim_match.range.end.try_into().unwrap(), After(false)),
            ))),
        };

        gb_data.features.push(gb_io::seq::Feature {
            kind: "primer_bind".into(),
            location,
            qualifiers: vec![("label".into(), name.to_owned().into())],
        });
    }

    let md = &data.metadata;

    gb_data.comments.clone_from(&md.comments);
    gb_data.source = Some(gb_io::seq::Source {
        source: md.source.clone().unwrap_or_default(),
        organism: md.organism.clone(),
    });

    if md.source.is_none() && md.organism.is_none() {
        gb_data.source = None;
    }

    // data.keywords = md.keywords.clone();
    gb_data.keywords = Some(md.plasmid_name.clone());
    gb_data.version.clone_from(&md.version);
    gb_data.accession.clone_from(&md.definition);
    gb_data.definition.clone_from(&md.definition);
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
