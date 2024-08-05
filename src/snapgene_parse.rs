//! Parse and write SnapGene DNA files. This converts between the Snapgene format for sequences, features,
//! primers, notes etc., and our own equivalents.
//!
//! (Unofficial file format description)[https://incenp.org/dvlpt/docs/binary-sequence-formats/binary-sequence-formats.pdf]

use std::{
    fs::{File, OpenOptions},
    io::{self, ErrorKind, Read, Seek, Write},
    path::Path,
    str,
};

use num_enum::TryFromPrimitive;
use quick_xml::de::from_str;
// use serde_xml_rs::{from_str, to_string};
use quick_xml::se::to_string;
use serde::Serialize;

use crate::{
    primer::{Primer, PrimerData},
    sequence::{
        seq_from_str, seq_to_str, Feature, FeatureDirection, FeatureType, Nucleotide, Seq,
        SeqTopology,
    },
    snapgene_parse::feature_xml::{
        FeatureSnapGene, Features, PrimerSnapGene, Primers, Qualifier, Segment,
    },
    util::{color_from_hex, color_to_hex},
    State,
};

const COOKIE_PACKET_LEN: usize = 14;

#[derive(Default)]
pub struct SnapgeneData {
    pub seq: Option<Seq>,
    pub topology: Option<SeqTopology>,
    pub features: Option<Vec<Feature>>,
    pub primers: Option<Vec<Primer>>,
}

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

/// DNA files are divided into packets. Packet structure:///
/// - A byte indicating the packet's type
/// - A BE 32-bit integer of packet len
/// - The payload
fn parse<R: Read + Seek>(file: &mut R) -> io::Result<(SnapgeneData)> {
    let buf = {
        let mut b = Vec::new();
        file.read_to_end(&mut b)?;
        b
    };

    let mut result = SnapgeneData::default();

    let mut i = 0;
    // Loop through each packet.
    loop {
        if i + 6 >= buf.len() {
            break;
        }
        let packet_type = match PacketType::try_from_primitive(buf[i]) {
            Ok(t) => t,
            Err(_) => PacketType::Unknown,
        };
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

                match parse_dna(&payload) {
                    Ok(v) => {
                        result.seq = Some(v.0);
                        result.topology = Some(v.1);
                    }
                    Err(e) => eprintln!("Error parsing DNA packet: {:?}", e),
                }
            }
            PacketType::Primers => match parse_primers(&payload) {
                Ok(v) => result.primers = Some(v),
                Err(e) => eprintln!("Error parsing Primers packet: {:?}", e),
            },
            PacketType::Notes => match parse_notes(&payload) {
                Ok(v) => (),
                Err(e) => eprintln!("Error parsing Notes packet: {:?}", e),
            },
            PacketType::Features => match parse_features(&payload) {
                Ok(v) => result.features = Some(v),
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
    use serde::{Deserialize, Serialize};

    #[derive(Debug, Serialize, Deserialize)]
    pub struct Features {
        #[serde(rename = "Feature", default)]
        pub inner: Vec<FeatureSnapGene>,
    }

    #[derive(Debug, Serialize, Deserialize)]
    pub struct FeatureSnapGene {
        #[serde(rename = "type", default)]
        pub feature_type: Option<String>,
        #[serde(rename = "Segment", default)]
        pub segments: Vec<Segment>,
        #[serde(rename = "Q", default)]
        pub qualifiers: Vec<Qualifier>,
        // pub qualifiers: Option<Vec<Qualifier>>,
        #[serde(default)]
        pub name: Option<String>,
        #[serde(default)]
        // pub directionality: Option<u8>,
        pub directionality: u8,
    }

    #[derive(Debug, Serialize, Deserialize)]
    pub struct Segment {
        #[serde(rename = "type", default)]
        pub segment_type: Option<String>,
        #[serde(default)]
        pub range: Option<String>,
        #[serde(default)]
        pub name: Option<String>,
        #[serde(default)]
        pub color: Option<String>, // Hex.
    }

    #[derive(Debug, Serialize, Deserialize)]
    pub struct Qualifier {
        #[serde(default)]
        pub name: String,
        #[serde(rename = "V", default)]
        pub values: Vec<QualifierValue>,
    }

    #[derive(Debug, Serialize, Deserialize)]
    pub struct QualifierValue {
        #[serde(default)]
        pub text: Option<String>,
        #[serde(default)]
        pub predef: Option<String>,
        #[serde(default)]
        pub int: Option<i32>,
    }

    #[derive(Debug, Serialize, Deserialize)]
    pub struct Primers {
        #[serde(rename = "Primer", default)]
        pub inner: Vec<PrimerSnapGene>,
    }

    // Note; We have left out the binding site and other fields, as they are not relevant for us at this time.
    #[derive(Debug, Serialize, Deserialize)]
    pub struct PrimerSnapGene {
        pub sequence: String,
        pub name: String,
        pub description: String, // Currently unused
    }

    #[derive(Debug, Serialize, Deserialize)]
    pub struct Notes {
        #[serde(rename = "Notes", default)]
        pub inner: Vec<Notes_>,
    }

    // Note; We have left out the binding site and other fields, as they are not relevant for us at this time.
    #[derive(Debug, Serialize, Deserialize)]
    pub struct Notes_ {}

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
            1 => FeatureDirection::Forward,
            2 => FeatureDirection::Reverse,
            // Some(1) => FeatureDirection::Forward,
            // Some(2) => FeatureDirection::Reverse,
            _ => FeatureDirection::None,
        };

        let feature_type = match &feature_sg.feature_type {
            Some(t) => FeatureType::from_external_str(t),
            None => FeatureType::default(),
        };

        for segment in &feature_sg.segments {
            let color_override = match &segment.color {
                Some(c) => color_from_hex(&c).ok(),
                None => None,
            };

            let index_range = match &segment.range {
                Some(r) => range_from_str(&r).unwrap_or((1, 1)),
                None => (1, 1),
            };

            result.push(Feature {
                index_range,
                feature_type,
                direction,
                label: name.clone(),
                color_override,
                notes: Default::default(), // todo: Add this.
            });
        }
        // todo: Handle qualifiers too?
    }

    Ok(result)
}

fn parse_primers(payload: &[u8]) -> io::Result<(Vec<Primer>)> {
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
        result.push(Primer {
            sequence: seq_from_str(&primer_sg.sequence),
            description: primer_sg.name.clone(),
        });
    }

    Ok(result)
}

fn parse_notes(payload: &[u8]) -> io::Result<(Vec<String>)> {
    let payload_str = str::from_utf8(payload).map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidData,
            format!("Unable to convert payload to string: {e}",),
        )
    })?;

    let notes: feature_xml::Notes = from_str(payload_str).map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidData,
            format!("Unable to parse notes: {e}"),
        )
    })?;

    let mut result = Vec::new();
    // for note in &notes.inner {
    // result.push(Primer {
    //     sequence: seq_from_str(&primer_sg.sequence),
    //     description: primer_sg.name.clone(),
    // });
    // }

    Ok(result)
}

fn range_from_str(range: &str) -> Result<(usize, usize), &'static str> {
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

    Ok((start, end))
}

/// Import data from a SnapGene .dna file into local state. This includes sequence, features, and primers.
pub fn import_snapgene(path: &Path) -> io::Result<(SnapgeneData)> {
    let mut file = File::open(path)?;
    parse(&mut file)
}

/// Add feature data to the buffer.
fn export_features(buf: &mut Vec<u8>, features: &[Feature]) -> io::Result<()> {
    let mut features_sg = Features { inner: Vec::new() };
    for feature in features {
        let directionality = match feature.direction {
            // FeatureDirection::Forward => Some(1),
            // FeatureDirection::Reverse => Some(2),
            // FeatureDirection::None => None,
            FeatureDirection::None => 0,
            FeatureDirection::Forward => 1,
            FeatureDirection::Reverse => 2,
        };

        let segments = vec![Segment {
            segment_type: None,
            range: Some(format!(
                "{}-{}",
                feature.index_range.0, feature.index_range.1
            )),
            name: None,
            color: feature.color_override.map(color_to_hex),
        }];

        features_sg.inner.push(FeatureSnapGene {
            feature_type: Some(feature.feature_type.to_string()),
            segments,
            qualifiers: Vec::new(),
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
fn export_primers(buf: &mut Vec<u8>, primers: &[PrimerData]) -> io::Result<()> {
    let mut primers_sg = Primers { inner: Vec::new() };
    for primer in primers {
        primers_sg.inner.push(PrimerSnapGene {
            sequence: seq_to_str(&primer.primer.sequence),
            name: primer.primer.description.clone(),
            description: "".to_owned(),
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
pub fn export_snapgene(
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

    let mut cookie_packet = [0; 19];

    cookie_packet[0] = PacketType::Cookie as u8;
    cookie_packet[1..5].clone_from_slice(&(COOKIE_PACKET_LEN as u32).to_be_bytes());
    cookie_packet[5..13].clone_from_slice(b"SnapGene");

    buf.extend(&cookie_packet);

    buf.push(PacketType::Dna as u8);
    buf.extend(((seq.len() + 1) as u32).to_be_bytes());

    let flag = match topology {
        SeqTopology::Circular => 1,
        SeqTopology::Linear => 0,
    };
    buf.push(flag);
    buf.extend(seq_to_str(seq).as_bytes());

    export_features(&mut buf, features)?;
    export_primers(&mut buf, primers)?;
    file.write_all(&buf)?;

    Ok(())
}
