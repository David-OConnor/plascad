//! Import and export GenBank files. This converts between this format for sequences, features,
//! primers, notes etc., and our own equivalents. The file format itself is handled by
//! [`bio_files::genbank`].
//!
//! (NIH article on GenBank)[https://www.ncbi.nlm.nih.gov/genbank/]

use std::{
    io::{self, ErrorKind},
    path::Path,
};

use bio_files::genbank::{
    GenBank, GenBankRecord, GenBankReference, KEY_ACCESSION, KEY_DATE, KEY_DIVISION, KEY_KEYWORDS,
    KEY_MOL_TYPE, KEY_ORGANISM, KEY_SOURCE, KEY_VERSION, format_date, parse_date,
};
use na_seq::{
    Nucleotide, SEQ_DESCRIPTION_KEY, SeqFeature, SeqRange, SeqTopology, Sequence, SequenceData,
    Strand, seq_complement,
};

use crate::{
    file_io::{GenericData, get_filename},
    misc_types::{Feature, FeatureDirection, FeatureType, Metadata, Reference},
    primer::{Primer, PrimerData, PrimerDirection, PrimerMatch},
    util::RangeIncl,
};

/// Read a file in the GenBank format. If it contains several records, this returns the first.
pub fn import_genbank(path: &Path) -> io::Result<GenericData> {
    let gb = GenBank::load(path)?;

    if gb.skipped_residues > 0 {
        eprintln!(
            "Left out {} residues from the GenBank file that can't be represented, e.g. N.",
            gb.skipped_residues
        );
    }

    let Some(record) = gb.records.into_iter().next() else {
        return Err(io::Error::new(
            ErrorKind::InvalidData,
            "No GenBank sequences found",
        ));
    };

    let GenBankRecord {
        seq: record_seq,
        references,
        comments,
    } = record;

    let seq = match record_seq.data {
        SequenceData::Dna(s) | SequenceData::Rna(s) => s,
        SequenceData::AminoAcid(_) => {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "This GenBank file contains a protein; only DNA and RNA are supported.",
            ));
        }
    };

    let topology = record_seq.topology.unwrap_or(SeqTopology::Linear);
    let (features, primers) = parse_features_primers(&record_seq.features, &seq);

    let meta = |key: &str| record_seq.metadata.get(key).cloned();

    let metadata = Metadata {
        plasmid_name: get_filename(path),
        date: meta(KEY_DATE).and_then(|d| parse_date(&d)),
        comments,
        definition: meta(SEQ_DESCRIPTION_KEY),
        molecule_type: meta(KEY_MOL_TYPE),
        division: meta(KEY_DIVISION).unwrap_or_default(),
        accession: meta(KEY_ACCESSION),
        version: meta(KEY_VERSION),
        keywords: meta(KEY_KEYWORDS),
        locus: record_seq.name.clone(),
        source: meta(KEY_SOURCE),
        organism: meta(KEY_ORGANISM),
        references: references.into_iter().map(reference_from_gb).collect(),
    };

    Ok(GenericData {
        seq,
        topology,
        features,
        primers,
        metadata,
    })
}

fn reference_from_gb(r: GenBankReference) -> Reference {
    Reference {
        description: r.description,
        authors: r.authors,
        consortium: r.consortium,
        title: r.title,
        journal: r.journal,
        pubmed: r.pubmed,
        remark: r.remark,
    }
}

fn reference_to_gb(r: &Reference) -> GenBankReference {
    GenBankReference {
        description: r.description.clone(),
        authors: r.authors.clone(),
        consortium: r.consortium.clone(),
        title: r.title.clone(),
        journal: r.journal.clone(),
        pubmed: r.pubmed.clone(),
        remark: r.remark.clone(),
    }
}

/// Parse features and primers, from GenBank's feature list.
fn parse_features_primers(
    features: &[SeqFeature],
    seq: &[Nucleotide],
) -> (Vec<Feature>, Vec<Primer>) {
    let mut result_ft = Vec::new();
    let mut primers = Vec::new();

    let compl = seq_complement(seq);

    for feature in features {
        let feature_type = FeatureType::from_external_str(&feature.kind);

        // The strand comes from the location being a complement, or the LEFT/RIGHT `direction`
        // qualifier. Primers and coding regions are always directional.
        let mut direction = match feature.strand {
            Strand::Forward => FeatureDirection::Forward,
            Strand::Reverse => FeatureDirection::Reverse,
            Strand::None => FeatureDirection::None,
        };

        match feature_type {
            FeatureType::Primer | FeatureType::CodingRegion => {
                if direction == FeatureDirection::None {
                    direction = FeatureDirection::Forward;
                }
            }
            _ => (),
        }

        let label = feature.label.clone();
        let notes = feature.qualifiers.clone();

        // We map multiple ranges, eg in the case of a GenBank `join` range type, to multiple features,
        // as our features currently only support a single range.
        for range in &feature.ranges {
            let mut range = RangeIncl::new(range.start, range.end);

            // Passing through the origin, most likely. Loop around past the end.
            if range.end < range.start {
                range = RangeIncl::new(range.start, seq.len() + range.end);
            }

            // GenBank stores primer bind sites (Which we treat as volatile), vice primer sequences.
            // Infer the sequence using the bind indices, and the main sequence.
            if feature_type == FeatureType::Primer {
                let sequence = match direction {
                    FeatureDirection::Reverse => {
                        let range = RangeIncl::new(
                            seq.len() - (range.end - 1),
                            seq.len() - (range.start - 1),
                        );

                        range.index_seq(&compl).unwrap_or_default()
                    }
                    _ => range.index_seq(seq).unwrap_or_default(),
                }
                .to_vec();

                // Perhaps improper way of storing primer descriptions
                let description = notes.first().map(|n| n.1.clone());

                let volatile = PrimerData::new(&sequence);
                primers.push(Primer {
                    sequence,
                    name: label.clone(),
                    description,
                    volatile,
                });
                continue;
            }

            result_ft.push(Feature {
                range,
                feature_type,
                direction,
                label: label.clone(),
                color_override: None,
                notes: notes.clone(),
            })
        }
    }

    (result_ft, primers)
}

fn strand_from_direction(direction: FeatureDirection) -> Strand {
    match direction {
        FeatureDirection::Forward => Strand::Forward,
        FeatureDirection::Reverse => Strand::Reverse,
        FeatureDirection::None => Strand::None,
    }
}

/// Export our local state into the GenBank format. This includes sequence, features, and primers.
pub fn export_genbank(
    data: &GenericData,
    primer_matches: &[(PrimerMatch, String)],
    path: &Path,
) -> io::Result<()> {
    let md = &data.metadata;

    let name = if md.locus.is_empty() {
        md.plasmid_name.clone()
    } else {
        md.locus.clone()
    };

    let mut seq = Sequence::new(SequenceData::Dna(data.seq.clone()), name);
    seq.topology = Some(data.topology);

    for feature in &data.features {
        seq.features.push(SeqFeature {
            kind: feature.feature_type.to_external_str(),
            label: feature.label.clone(),
            ranges: vec![SeqRange::new(feature.range.start, feature.range.end)],
            strand: strand_from_direction(feature.direction),
            color: None,
            qualifiers: feature.notes.clone(),
        });
    }

    for (prim_match, name) in primer_matches {
        let mut qualifiers = Vec::new();

        // todo: This is a sloppy way of accessing the primer.
        if let Some(descrip) = data
            .primers
            .iter()
            .find(|p| p.name == *name)
            .and_then(|p| p.description.as_ref())
        {
            qualifiers.push(("note".to_owned(), descrip.to_owned()));
        }

        // todo: Make sure we're not getting a duplicate label.
        seq.features.push(SeqFeature {
            kind: "primer_bind".to_owned(),
            label: name.clone(),
            ranges: vec![SeqRange::new(prim_match.range.start, prim_match.range.end)],
            strand: match prim_match.direction {
                PrimerDirection::Forward => Strand::Forward,
                PrimerDirection::Reverse => Strand::Reverse,
            },
            color: None,
            qualifiers,
        });
    }

    let mut meta = |key: &str, val: Option<&String>| {
        if let Some(v) = val.filter(|v| !v.is_empty()) {
            seq.metadata.insert(key.to_owned(), v.clone());
        }
    };

    meta(SEQ_DESCRIPTION_KEY, md.definition.as_ref());
    meta(KEY_MOL_TYPE, md.molecule_type.as_ref());
    meta(KEY_DIVISION, Some(&md.division));
    meta(KEY_ACCESSION, md.accession.as_ref());
    meta(KEY_VERSION, md.version.as_ref());
    meta(KEY_KEYWORDS, md.keywords.as_ref());
    meta(KEY_SOURCE, md.source.as_ref());
    meta(KEY_ORGANISM, md.organism.as_ref());

    let date = md.date.map(|(y, m, d)| format_date(y, m, d));
    meta(KEY_DATE, date.as_ref());

    GenBank {
        records: vec![GenBankRecord {
            seq,
            references: md.references.iter().map(reference_to_gb).collect(),
            comments: md.comments.clone(),
        }],
        skipped_residues: 0,
    }
    .save(path)
}
