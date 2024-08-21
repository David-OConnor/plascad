//! Parse and write SnapGene DNA files. This converts between the Snapgene format for sequences, features,
//! primers, notes etc., and our own equivalents.
//!
//! (Unofficial file format description)[https://incenp.org/dvlpt/docs/binary-sequence-formats/binary-sequence-formats.pdf]
//! DNA files are divided into packets. Packet structure:
//! - A byte indicating the packet's type
//! - A big endian 32-bit unsigned integer of packet len
//! - The payload

use std::{
    collections::HashMap,
    fs::File,
    io::{self, ErrorKind, Read, Write},
    ops::RangeInclusive,
    path::Path,
    str,
};

use num_enum::TryFromPrimitive;
use quick_xml::{de::from_str, se::to_string};

// We remove these from SnapGene feature qualifiers.
const HTML_TAGS: [&str; 8] = [
    "<html>", "</html>", "<body>", "</body>", "<i>", "</i>", "<b>", "</b>",
];

use crate::{
    file_io::{
        get_filename,
        snapgene::feature_xml::{
            FeatureSnapGene, Features, Notes, PrimerSnapGene, Primers, Qualifier, QualifierValue,
            Segment,
        },
        GenericData,
    },
    primer::{Primer, PrimerData},
    sequence::{
        seq_from_str, seq_to_str, Feature, FeatureDirection, FeatureType, Nucleotide, Seq,
        SeqTopology,
    },
    util::{color_from_hex, color_to_hex, RangeIncl},
};

const COOKIE_PACKET_LEN: usize = 14;

#[derive(Debug, TryFromPrimitive)]
#[repr(u8)]
// todo: We are observing other packet types: 0x9, 0x3, 0x11, 0x8, 0xd, 0xe, 0x1c
enum PacketType {
    /// The cookie is always the first packet.
    Cookie = 0x09,
    Dna = 0x00,
    Primers = 0x05,
    Notes = 0x06,
    Features = 0x0a,
    AdditionalSequenceProperties = 0x08,
    AlignableSequences = 0x11,
    CustomEnzymeSets = 0xe,
    Unknown = 0x99, // Placeholder for encountering one we don't recognize
}

/// Import a file in SnapGene's DNA format into local state. This includes sequence, features, and primers.
pub fn import_snapgene(path: &Path) -> io::Result<GenericData> {
    let mut file = File::open(path)?;

    let buf = {
        let mut b = Vec::new();
        file.read_to_end(&mut b)?;
        b
    };

    let mut result = GenericData::default();

    result.metadata.plasmid_name = get_filename(path);

    let mut i = 0;

    loop {
        if i + 6 >= buf.len() {
            break;
        }
        let packet_type = PacketType::try_from_primitive(buf[i]).unwrap_or(PacketType::Unknown);
        i += 1;

        let payload_len = u32::from_be_bytes(buf[i..i + 4].try_into().unwrap()) as usize;
        i += 4;

        if i + payload_len + 1 > buf.len() {
            eprintln!(
                "Error parsing DNA file: Payload would exceed file length. Index: {}, buf len: {}",
                i + payload_len,
                buf.len()
            );
            break;
        }

        let payload = &buf[i..i + payload_len];
        i += payload_len;

        match packet_type {
            PacketType::Cookie => {
                if payload_len != COOKIE_PACKET_LEN {
                    eprintln!("Invalid cookie packet length: {}", payload_len);
                }
                if &payload[..8] != b"SnapGene" {
                    eprintln!("Invalid cookie payload: {:?}", &payload[0..8]);
                }
                // The next bytes describe the type of seq (1 for DNA, export version, and import version)
                // at indexes 0xd, 0xf, and 0x11 respectively.
            }
            PacketType::Dna => {
                // todo: Note: This doesn't properly handle if there are multiple DNA packets.
                // todo: How should we do that?

                match parse_dna(payload) {
                    Ok(v) => {
                        result.seq = v.0;
                        result.topology = v.1;
                    }
                    Err(e) => eprintln!("Error parsing DNA packet: {:?}", e),
                }
            }
            PacketType::Primers => match parse_primers(payload) {
                Ok(v) => result.primers = v,
                Err(e) => eprintln!("Error parsing Primers packet: {:?}", e),
            },
            PacketType::Notes => match parse_notes(payload) {
                Ok(_v) => {
                    // if !v.inner.is_empty() {
                    //     // todo: Are there ever multiple notes?
                    //     result.metadata.plasmid_name = v.inner[0].title.clone();
                    // }
                    // todo: Other fields, references etc. Compare in SnapGene itself, and how snapgene
                    // todo parses and exports a GenBank file.
                }
                Err(e) => eprintln!("Error parsing Notes packet: {:?}", e),
            },
            PacketType::Features => match parse_features(payload) {
                Ok(v) => result.features = v,
                Err(e) => eprintln!("Error parsing Features packet: {:?}", e),
            },
            PacketType::Unknown => {
                println!(
                    "Unknown packet type: {:x} len: {payload_len}",
                    buf[i - 5 - payload_len]
                );

                let payload_str = str::from_utf8(payload)
                    .map_err(|e| {
                        io::Error::new(
                            ErrorKind::InvalidData,
                            format!("Unable to convert payload to string: {e}",),
                        )
                    })
                    .ok();

                println!("Payload str: \n{:?}", payload_str);
            }
            _ => (),
        }
    }

    Ok(result)
}

fn parse_dna(payload: &[u8]) -> io::Result<(Seq, SeqTopology)> {
    if payload.is_empty() {
        return Err(io::Error::new(ErrorKind::InvalidData, "Empty DNA packet"));
    }

    let flags = payload[0];
    let sequence = &payload[1..];

    let mut seq = Vec::new();

    for nt in sequence {
        match Nucleotide::from_u8_letter(*nt) {
            Ok(n) => seq.push(n),
            Err(_) => {
                eprintln!("Unexpected char in DNA sequence: {:?}", nt);
            }
        }
    }

    let topology = if flags & 0x01 != 0 {
        SeqTopology::Circular
    } else {
        SeqTopology::Linear
    };

    println!("Flags: {flags}");

    Ok((seq, topology))
}

// todo: Consider a sub-module for XML parsing.
mod feature_xml {
    use std::str::FromStr;

    use serde::{Deserialize, Deserializer, Serialize};

    #[derive(Debug, Serialize, Deserialize)]
    pub struct Features {
        #[serde(rename = "Feature", default)]
        pub inner: Vec<FeatureSnapGene>,
    }

    // Workaround for parsing "" into None not supported natively by quick-xml/serde.
    fn deserialize_directionality<'de, D>(deserializer: D) -> Result<Option<u8>, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s: Option<String> = Option::deserialize(deserializer)?;
        if let Some(s) = s {
            if s.is_empty() {
                Ok(None)
            } else {
                u8::from_str(&s).map(Some).map_err(serde::de::Error::custom)
            }
        } else {
            Ok(None)
        }
    }

    #[derive(Debug, Serialize, Deserialize)]
    pub struct FeatureSnapGene {
        #[serde(rename = "@type")]
        pub feature_type: Option<String>,
        #[serde(
            rename = "@directionality",
            deserialize_with = "deserialize_directionality",
            default
        )]
        pub directionality: Option<u8>,
        #[serde(rename = "@name", default)]
        pub name: Option<String>,
        // Other Feature attributes: allowSegmentOverlaps (0/1), consecutiveTranslationNumbering (0/1)
        #[serde(rename = "Segment", default)]
        pub segments: Vec<Segment>,
        #[serde(rename = "Q", default)]
        pub qualifiers: Vec<Qualifier>,
    }

    #[derive(Debug, Serialize, Deserialize)]
    pub struct Segment {
        #[serde(rename = "@type", default)]
        pub segment_type: Option<String>,
        #[serde(rename = "@range", default)]
        pub range: Option<String>,
        #[serde(rename = "@name", default)]
        pub name: Option<String>,
        #[serde(rename = "@color", default, skip_serializing_if = "Option::is_none")]
        pub color: Option<String>, // Hex.
                                   // Other fields: "translated": 0/1
    }

    #[derive(Debug, Serialize, Deserialize)]
    pub struct Qualifier {
        #[serde(rename = "@name")]
        pub name: String,
        #[serde(rename = "V", default)]
        pub values: Vec<QualifierValue>,
    }

    #[derive(Debug, Serialize, Deserialize)]
    pub struct QualifierValue {
        #[serde(rename = "@text", default, skip_serializing_if = "Option::is_none")]
        pub text: Option<String>,
        #[serde(rename = "@predef", default, skip_serializing_if = "Option::is_none")]
        pub predef: Option<String>,
        #[serde(rename = "@int", default, skip_serializing_if = "Option::is_none")]
        pub int: Option<i32>,
    }

    #[derive(Debug, Serialize, Deserialize)]
    pub struct Primers {
        #[serde(rename = "Primer", default)]
        pub inner: Vec<PrimerSnapGene>,
    }

    // Note; We have left out the binding site and a number of other fields, as they are not relevant
    // for us at this time. This also includes melting temperature, which we calculate.
    #[derive(Debug, Serialize, Deserialize)]
    pub struct PrimerSnapGene {
        #[serde(rename = "@sequence")]
        pub sequence: String,
        #[serde(rename = "@name")]
        pub name: String,
        #[serde(rename = "@description")]
        pub description: String,
    }

    #[derive(Debug, Serialize, Deserialize)]
    pub struct Notes {
        #[serde(rename = "Notes", default)]
        pub inner: Vec<Notes_>,
    }

    // Note; We have left out the binding site and other fields, as they are not relevant for us at this time.
    #[derive(Debug, Serialize, Deserialize)]
    pub struct Notes_ {
        #[serde(rename = "UUID")]
        pub uuid: String,
        #[serde(rename = "Type")]
        pub type_: String,
        #[serde(rename = "ConfirmedExperimentally")]
        pub confirmed_experimentally: u8, // todo? `0` in example
        #[serde(rename = "CreatedBy")]
        pub created_by: String,
        #[serde(rename = "SequenceClass")]
        pub sequence_class: String,
        #[serde(rename = "TransformedInto")]
        pub transformed_into: String,
        // todo: How do we handle LastModified and Created, given they have timestamps in the field name?
        #[serde(rename = "Reference journal")]
        pub reference_journal: String,
        pub doi: String,
        pub pages: String,
        #[serde(rename = "pubMedID")]
        pub pub_med_id: String,
        pub title: String,
        pub date: String,
        pub authors: String,
        #[serde(rename = "journalName")]
        pub nournal_name: String,
        pub volume: String,
    }

    // <Notes>
    //     <UUID>0962493c-08f0-4964-91b9-24840fea051e</UUID>
    //     <Type>Synthetic</Type>
    //     <ConfirmedExperimentally>0</ConfirmedExperimentally>
    //     <Created UTC="22:48:15">2024.7.12</Created>
    //     <LastModified UTC="18:2:52">2024.7.14</LastModified>
    //     <CreatedBy>SnapGene License</CreatedBy>
    //     <SequenceClass>UNA</SequenceClass>
    //     <TransformedInto>DH5α™</TransformedInto>
    // </Notes>
}

fn parse_features(payload: &[u8]) -> io::Result<Vec<Feature>> {
    let payload_str = str::from_utf8(payload).map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidData,
            format!("Unable to convert payload to string: {e}",),
        )
    })?;

    // println!("\n\n\nPayload str: {:?}\n\n\n", payload_str);

    let features: Features = from_str(payload_str).map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidData,
            format!("Unable to parse features: {e}"),
        )
    })?;

    let mut result = Vec::new();

    for feature_sg in &features.inner {
        // Note: Our model does not include the concept of segments; treat each SnapGene segment as a new feature.
        let name = feature_sg.name.clone().unwrap_or(String::new());

        let direction = match feature_sg.directionality {
            Some(1) => FeatureDirection::Forward,
            Some(2) => FeatureDirection::Reverse,
            _ => FeatureDirection::None,
        };

        let feature_type = match &feature_sg.feature_type {
            Some(t) => FeatureType::from_external_str(t),
            None => FeatureType::default(),
        };

        // let mut notes = HashMap::new();
        let mut notes = Vec::new();
        for qual in &feature_sg.qualifiers {
            // It seems there is generally one value per qualifier. In the case there are multiple,
            // we will parse them as separate notes.
            for val in &qual.values {
                // Generally, each val only has one of text, int, predef
                let mut v = String::new();
                if let Some(t) = &val.int {
                    v = t.to_string();
                }

                if let Some(t) = &val.predef {
                    v.clone_from(t);
                }
                if let Some(t) = &val.text {
                    v.clone_from(t);
                }

                // Remove the HTML tags and related that SnapGene inserts into qual values.
                for html_tag in HTML_TAGS {
                    v = v.replace(html_tag, "");
                }

                // notes.insert(qual.name.clone(), v);
                notes.push((qual.name.clone(), v));
            }
        }

        // Note: We currently parse multiple segments as separate features, but SnapGene has the concept
        // of multiple segments per feature. These share all data except the <Segment tag attributes.
        // (name, range, color, type, translated etc)
        for segment in &feature_sg.segments {
            let color_override = match &segment.color {
                Some(c) => color_from_hex(c).ok(),
                None => None,
            };

            let range = match &segment.range {
                Some(r) => range_from_str(r).unwrap_or(RangeIncl::new(1, 1)),
                None => RangeIncl::new(1, 1),
            };

            result.push(Feature {
                range,
                feature_type,
                direction,
                label: name.clone(),
                color_override,
                notes: notes.clone(),
            });
        }
    }

    Ok(result)
}

fn parse_primers(payload: &[u8]) -> io::Result<Vec<Primer>> {
    let payload_str = str::from_utf8(payload).map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidData,
            format!("Unable to convert payload to string: {e}",),
        )
    })?;

    let primers: Primers = from_str(payload_str).map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidData,
            format!("Unable to parse primers: {e}"),
        )
    })?;

    let mut result = Vec::new();
    for primer_sg in &primers.inner {
        let seq = seq_from_str(&primer_sg.sequence);
        let volatile = PrimerData::new(&seq);
        result.push(Primer {
            sequence: seq,
            name: primer_sg.name.clone(),
            description: Some(primer_sg.description.clone()),
            volatile,
        });
    }

    Ok(result)
}

// fn parse_notes(payload: &[u8]) -> io::Result<Vec<String>> {
fn parse_notes(payload: &[u8]) -> io::Result<Notes> {
    let payload_str = str::from_utf8(payload).map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidData,
            format!("Unable to convert payload to string: {e}",),
        )
    })?;

    println!("Notes string: \n\n{:?}\n\n", payload_str);

    // todo: Is this a strict format, or arbitary notes?

    let notes: Notes = from_str(payload_str).map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidData,
            format!("Unable to parse notes: {e}"),
        )
    })?;
    let result = notes;

    Ok(result)
}

fn range_from_str(range: &str) -> Result<RangeIncl, &'static str> {
    let parts: Vec<&str> = range.split('-').collect();
    if parts.len() != 2 {
        return Err("Invalid range format");
    }

    let start = parts[0]
        .parse::<usize>()
        .map_err(|_| "Invalid number in range")?;
    let end = parts[1]
        .parse::<usize>()
        .map_err(|_| "Invalid number in range")?;

    Ok(RangeIncl::new(start, end))
}

/// Add feature data to the buffer.
fn export_features(buf: &mut Vec<u8>, features: &[Feature]) -> io::Result<()> {
    let mut features_sg = Features { inner: Vec::new() };
    for feature in features {
        let directionality = match feature.direction {
            FeatureDirection::Forward => Some(1),
            FeatureDirection::Reverse => Some(2),
            FeatureDirection::None => None,
        };

        let segments = vec![Segment {
            segment_type: None,
            range: Some(format!("{}-{}", feature.range.start, feature.range.end)),
            name: None,
            color: feature.color_override.map(color_to_hex),
        }];

        let mut qualifiers = Vec::new();
        for (key, value) in &feature.notes {
            // todo: If int parsable, consider saving as an int.
            qualifiers.push(Qualifier {
                name: key.to_string(),
                values: vec![QualifierValue {
                    text: Some(value.to_string()),
                    predef: None,
                    int: None,
                }],
            })
        }

        features_sg.inner.push(FeatureSnapGene {
            feature_type: Some(feature.feature_type.to_string()),
            segments,
            qualifiers,
            name: Some(feature.label.clone()),
            directionality,
        });
    }

    let xml_str = to_string(&features_sg).map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidData,
            format!("Unable to convert features to an XML string: {e}"),
        )
    })?;

    let xml = xml_str.into_bytes();

    buf.push(PacketType::Features as u8);
    buf.extend((xml.len() as u32).to_be_bytes());
    buf.extend(&xml);

    Ok(())
}

/// Add primer data to the buffer.
fn export_primers(buf: &mut Vec<u8>, primers: &[Primer]) -> io::Result<()> {
    let mut primers_sg = Primers { inner: Vec::new() };
    for primer in primers {
        primers_sg.inner.push(PrimerSnapGene {
            sequence: seq_to_str(&primer.sequence),
            name: primer.name.clone(),
            description: primer.description.clone().unwrap_or_default(),
        });
    }

    let xml_str = to_string(&primers_sg).map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidData,
            format!("Unable to convert primers to an XML string: {e}"),
        )
    })?;

    let xml = xml_str.into_bytes();

    buf.push(PacketType::Features as u8);
    buf.extend((xml.len() as u32).to_be_bytes());
    buf.extend(&xml);

    Ok(())
}

/// Export our local state into the SnapGene dna format. This includes sequence, features, and primers.
pub fn export_snapgene(data: &GenericData, path: &Path) -> io::Result<()> {
    let mut file = File::create(path)?;

    let mut buf = Vec::new();

    let mut cookie_packet = [0; 19];

    cookie_packet[0] = PacketType::Cookie as u8;
    cookie_packet[1..5].clone_from_slice(&(COOKIE_PACKET_LEN as u32).to_be_bytes());
    cookie_packet[5..13].clone_from_slice(b"SnapGene");

    buf.extend(&cookie_packet);

    buf.push(PacketType::Dna as u8);
    buf.extend(((data.seq.len() + 1) as u32).to_be_bytes());

    let flag = match data.topology {
        SeqTopology::Circular => 1,
        SeqTopology::Linear => 0,
    };
    buf.push(flag);
    buf.extend(seq_to_str(&data.seq).as_bytes());

    export_features(&mut buf, &data.features)?;
    export_primers(&mut buf, &data.primers)?;
    file.write_all(&buf)?;

    Ok(())
}
