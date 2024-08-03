//! Parse and write SnapGene DNA files. This converts between the Snapgene format for sequences, features,
//! primers, notes etc, and our own equivalents.
//!
//! (Unofficial file format description)[https://incenp.org/dvlpt/docs/binary-sequence-formats/binary-sequence-formats.pdf]

use std::{
    fs::File,
    io::{self, ErrorKind, Read, Seek},
    path::Path,
    str,
};
use std::collections::HashMap;
use serde_xml_rs::from_str;
// use xml::reader::{EventReader, XmlEvent};

use chrono::NaiveDate;

use crate::{
    primer::{Primer, PrimerData},
    sequence::{Feature, Nucleotide, Seq, SeqTopology},
};
use crate::sequence::{FeatureDirection, FeatureType};
use crate::util::color_from_hex;

#[derive(Default)]
pub struct SnapgeneData {
    pub seq: Option<Seq>,
    pub topology: Option<SeqTopology>,
    pub features: Option<Vec<Feature>>,
    pub primers: Option<Vec<Primer>>,
}
//
// // #[derive(Debug, Deserialize)]
// struct SeqRecord {
//     pub seq: Option<String>,
//     pub annotations: Annotations,
//     pub id: Option<String>,
//     pub name: Option<String>,
//     pub  description: Option<String>,
//     pub features: Vec<SeqFeature>,
// }
//
// // #[derive(Debug, Deserialize)]
// struct Annotations {
//     pub molecule_type: Option<String>,
//     pub topology: Option<String>,
//     pub data_file_division: Option<String>,
//     pub date: Option<NaiveDate>,
// }
//
// // #[derive(Debug, Deserialize)]
// struct SeqFeature {
//     pub location: SimpleLocation,
//     pub feature_type: String,
//     pub qualifiers: Qualifiers,
// }
//
//
// // #[derive(Debug, Deserialize)]
// struct Qualifiers {
//     pub parts: Option<String>,
//     pub label: Option<Vec<String>>,
//     pub name: Option<Vec<String>>,
// }

#[derive(Debug)]
#[repr(u8)]
// todo: We are observing other packet types: 0x9, 0x3, 0x11, 0x8, 0xd, 0xe, 0x1c
enum PacketType {
    /// The cookie is always the first packet.
    Cookie = 0x09,
    Dna = 0x00,
    Primers = 0x05,
    Notes = 0x06,
    Features = 0x0A,
    Unknown = 0x99, // Placeholder for encountering one we don't recognize
}

impl PacketType {
    pub fn from_byte(byte: u8) -> io::Result<Self> {
        match byte {
            0x09 => Ok(Self::Cookie),
            0x00 => Ok(Self::Dna),
            0x05 => Ok(Self::Primers),
            0x06 => Ok(Self::Notes),
            0x0A => Ok(Self::Features),
            _ => Err(io::Error::new(
                ErrorKind::InvalidData,
                "Invalid byte for packet SnapGene packet type",
            )),
        }
    }
}

fn read_bytes<R: Read>(reader: &mut R, length: usize) -> io::Result<Vec<u8>> {
    let mut buffer = vec![0; length];
    reader.read_exact(&mut buffer)?;
    Ok(buffer)
}

/// DNA files are divided into packets. Packet structure:///
/// - A byte indicating the packet's type
/// - A BE 32-bit integer of packet len
/// - The payload
fn parse<R: Read + Seek>(file: &mut R) -> io::Result<(SnapgeneData)> {
    println!("Starting read...");

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
        let packet_type = match PacketType::from_byte(buf[i]) {
            Ok(t) => t,
            Err(_) => {
                println!("Unknown packet type: {:x}", buf[i]);
                PacketType::Unknown
            }
        };
        i += 1;

        let payload_len = u32::from_be_bytes(buf[i..i + 4].try_into().unwrap()) as usize;
        i += 4;

        println!("\nPacket type: {:?}, len: {:?}", packet_type, payload_len);

        if i + payload_len + 1 > buf.len() {
            eprintln!(
                "Error parsing DNA file: Payload would exceed file length. Index: {}, buff len: {}",
                i + payload_len,
                buf.len()
            );
            break;
        }

        // let end = min(payload_len + i, buf.len() - 1); // todo temp. want ot see the dna eq
        let payload = &buf[i..i + payload_len];
        // let payload = &buf[i..end];
        i += payload_len;

        match packet_type {
            PacketType::Cookie => {
                println!("Cookie packet type.");
                if payload_len != 14 {
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

                println!("DNA packet type");
                match parse_dna_packet(&payload) {
                    Ok(v) => {
                        result.seq = Some(v.0);
                        result.topology = Some(v.1);
                    }
                    Err(e) => eprintln!("Error parsing DNA packet: {:?}", e),
                }
            }
            PacketType::Primers => {
                println!("Primers packet type");
            }
            PacketType::Notes => {
                println!("Notes packet type");
            }
            PacketType::Features => {
                println!("Features packet type");
                match parse_features_packet(&payload) {
                    Ok(v) => result.features = Some(v),
                    Err(e) => eprintln!("Error parsing Features packet: {:?}", e),
                }
            }
            PacketType::Unknown => {
                println!("Unknown packet type: {:?}", buf[i - 5 - payload_len]);
            }
        }
    }

    Ok(result)
    //
    //
    // let length = read_bytes(buffer, 4)?;
    // let length = u32::from_be_bytes(length.try_into().unwrap()) as usize;
    //
    // let payload = read_bytes(buffer, length)?;
}

fn parse_dna_packet(payload: &[u8]) -> io::Result<(Seq, SeqTopology)> {
    if payload.is_empty() {
        return Err(io::Error::new(ErrorKind::InvalidData, "Empty DNA packet"));
    }

    let flags = payload[0];
    let sequence = &payload[1..];

    let mut seq = Vec::new();

    for nt in sequence {
        match Nucleotide::from_u8(*nt) {
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

fn parse_primers_packet(payload: &[u8]) -> io::Result<()> {
    Ok(())
}

fn parse_notes_packet(payload: &[u8]) -> io::Result<()> {
    let xml = str::from_utf8(payload).unwrap();
    // let notes: Notes = from_str(xml).unwrap();

    // Process notes data
    Ok(())
}

// // todo: Temp
// struct Record {
//     features: Vec<SeqFeature>,
// }

struct SimpleLocation {
    start: usize,
    end: usize,
    strand: char,
}

impl SimpleLocation {
    fn new(start: usize, end: usize, strand: char) -> Self {
        SimpleLocation { start, end, strand }
    }

    fn add(&self, other: SimpleLocation) -> SimpleLocation {
        // Combine two SimpleLocations, assuming strand and other properties are managed appropriately
        // This is a placeholder, actual implementation may differ
        SimpleLocation::new(self.start, other.end, self.strand)
    }
}

fn parse_location(rangespec: &str, strand: char, record_len: usize, is_primer: bool) -> SimpleLocation {
    let mut parts = rangespec.split('-');
    let start: usize = parts.next().unwrap().parse().unwrap();
    let end: usize = parts.next().unwrap().parse().unwrap();

    // Account for SnapGene's 1-based coordinates
    let mut start = start - 1;
    let mut end = end;

    if is_primer {
        // Primers' coordinates in SnapGene files are shifted by -1 for some reasons
        start += 1;
        end += 1;
    }

    if start > end {
        // Range wrapping the end of the sequence
        let l1 = SimpleLocation::new(start, record_len, strand);
        let l2 = SimpleLocation::new(0, end, strand);
        l1.add(l2)
    } else {
        SimpleLocation::new(start, end, strand)
    }
}

// todo: Consider a sub-module for XML parsing.
mod feature_xml {
    use std::collections::HashMap;
    use serde::Deserialize;

    #[derive(Debug, Deserialize)]
    pub struct Features {
        #[serde(rename = "Feature", default)]
        pub features: Vec<FeatureSnapGene>,
    }

    #[derive(Debug, Deserialize)]
    pub struct FeatureSnapGene {
        #[serde(rename = "type")]
        pub feature_type: Option<String>,
        #[serde(rename = "Segment", default)]
        pub segments: Vec<Segment>,
        #[serde(rename = "Q", default)]
        pub qualifiers: Vec<Qualifier>,
        pub name: Option<String>,
        pub directionality: Option<u8>,
    }

    #[derive(Debug, Deserialize)]
    pub struct Segment {
        #[serde(rename = "type")]
        pub segment_type: Option<String>,
        pub range: Option<String>,
        pub name: Option<String>,
        pub color: Option<String>, // Hex.
    }

    #[derive(Debug, Deserialize)]
    pub struct Qualifier {
        pub name: String,
        #[serde(rename = "V", default)]
        pub values: Vec<QualifierValue>,
    }

    #[derive(Debug, Deserialize)]
    pub struct QualifierValue {
        pub text: Option<String>,go +
        pub predef: Option<String>,
        pub int: Option<i32>,
    }

    #[derive(Debug)]
    pub struct SeqFeature {
        pub location: Option<String>,
        pub feature_type: String,
        pub qualifiers: HashMap<String, Vec<String>>,
    }
}

fn range_from_str(range: &str) -> Result<(usize, usize), &'static str> {
    let parts: Vec<&str> = range.split('-').collect();
    if parts.len() != 2 {
        return Err("Invalid range format");
    }

    let start = parts[0].parse::<usize>().map_err(|_| "Invalid number in range")?;
    let end = parts[1].parse::<usize>().map_err(|_| "Invalid number in range")?;

    Ok((start, end))
}

fn parse_features_packet(payload: &[u8]) -> io::Result<Vec<Feature>> {
    let payload_str = str::from_utf8(payload).unwrap(); // todo: Handle

    let features: feature_xml::Features = from_str(payload_str).unwrap(); // todo:  Handle

    let mut result = Vec::new();

    for feature_sg in &features.features {
        // Note: Our model does not include the concept of segments; treat each SnapGene segment as a new feature.
        let name = feature_sg.name.clone().unwrap_or(String::new());

        println!("feature name raw: {:?}", feature_sg.name);

        let direction = match feature_sg.directionality {
            Some(1) => FeatureDirection::Forward,
            Some(2) => FeatureDirection::Reverse,
            _ => FeatureDirection::None,
        };

        let feature_type = match &feature_sg.feature_type {
            Some(t) => FeatureType::from_snapgene_str(t),
            None => FeatureType::default()
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
            });
        }
    }

    Ok(result)
}

// #[derive(Debug, Deserialize)]
struct Notes {
    type_: Option<String>,
    last_modified: Option<String>,
    accession_number: Option<String>,
    comments: Option<String>,
}


pub fn import_snapgene(path: &Path) -> io::Result<(SnapgeneData)> {
    let mut file = File::open(path)?;
    parse(&mut file)
}
