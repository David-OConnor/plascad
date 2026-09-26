//! Import and export SnapGene DNA files. This converts between the SnapGene format for sequences,
//! features, primers, notes etc., and our own equivalents. The file format itself is handled by
//! [`bio_files::snapgene`].

use std::{
    io::{self, ErrorKind},
    path::Path,
};

use bio_files::snapgene::{SnapGene, SnapGenePrimer};
use na_seq::{SeqFeature, SeqRange, Sequence, SequenceData, Strand};

use crate::{
    file_io::{GenericData, get_filename},
    misc_types::{Feature, FeatureDirection, FeatureType},
    primer::{Primer, PrimerData},
    util::{RangeIncl, color_from_hex, color_to_hex},
};

/// Import a file in SnapGene's DNA format into local state. This includes sequence, features, and primers.
pub fn import_snapgene(path: &Path) -> io::Result<GenericData> {
    let sg = SnapGene::load(path)?;

    if sg.skipped_residues > 0 {
        eprintln!(
            "Left out {} residues from the SnapGene file that can't be represented, e.g. N.",
            sg.skipped_residues
        );
    }

    let seq = match &sg.seq.data {
        SequenceData::Dna(s) => s.clone(),
        _ => {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "The SnapGene file doesn't contain DNA",
            ));
        }
    };

    let mut result = GenericData {
        seq,
        topology: sg.seq.topology.unwrap_or_default(),
        ..Default::default()
    };

    result.metadata.plasmid_name = get_filename(path);
    result.metadata.definition = sg.seq.description().map(str::to_owned);

    for feature in &sg.seq.features {
        let direction = match feature.strand {
            Strand::Forward => FeatureDirection::Forward,
            Strand::Reverse => FeatureDirection::Reverse,
            Strand::None => FeatureDirection::None,
        };

        let color_override = feature
            .color
            .as_deref()
            .and_then(|c| color_from_hex(c).ok());

        // Our model does not include the concept of segments; treat each SnapGene segment as a
        // new feature.
        for range in &feature.ranges {
            result.features.push(Feature {
                range: RangeIncl::new(range.start, range.end),
                feature_type: FeatureType::from_external_str(&feature.kind),
                direction,
                label: feature.label.clone(),
                color_override,
                notes: feature.qualifiers.clone(),
            });
        }
    }

    result.primers = sg
        .primers
        .into_iter()
        .map(|p| {
            let volatile = PrimerData::new(&p.sequence);
            Primer {
                sequence: p.sequence,
                name: p.name,
                description: p.description,
                volatile,
            }
        })
        .collect();

    Ok(result)
}

/// Export our local state into the SnapGene dna format. This includes sequence, features, and primers.
pub fn export_snapgene(data: &GenericData, path: &Path) -> io::Result<()> {
    let mut seq = Sequence::new(
        SequenceData::Dna(data.seq.clone()),
        data.metadata.plasmid_name.clone(),
    );
    seq.topology = Some(data.topology);

    seq.features = data
        .features
        .iter()
        .map(|feature| SeqFeature {
            kind: feature.feature_type.to_external_str(),
            label: feature.label.clone(),
            ranges: vec![SeqRange::new(feature.range.start, feature.range.end)],
            strand: match feature.direction {
                FeatureDirection::Forward => Strand::Forward,
                FeatureDirection::Reverse => Strand::Reverse,
                FeatureDirection::None => Strand::None,
            },
            color: feature.color_override.map(color_to_hex),
            qualifiers: feature.notes.clone(),
        })
        .collect();

    let primers = data
        .primers
        .iter()
        .map(|p| SnapGenePrimer {
            sequence: p.sequence.clone(),
            name: p.name.clone(),
            description: p.description.clone(),
        })
        .collect();

    SnapGene {
        seq,
        primers,
        skipped_residues: 0,
    }
    .save(path)
}
