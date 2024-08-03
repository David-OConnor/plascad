//! Parse and write SnapGene DNA files. This converts between the Snapgene format for sequences, features,
//! primers, notes etc, and our own equivalents.

use std::{
    fs::File,
    io::{self, ErrorKind, Read, Seek, SeekFrom},
    path::Path,
    str,
};

use chrono::NaiveDate;

use crate::{
    primer::{Primer, PrimerData},
    sequence::{Feature, Nucleotide, Seq, SeqTopology},
};
// use serde::Deserialize;
// use serde_xml_rs::from_str;

#[derive(Default)]
pub struct SnapgeneData {
    pub seq: Option<Seq>,
    pub topology: Option<SeqTopology>,
    pub features: Option<Vec<Feature>>,
    pub primers: Option<Vec<Primer>>,
}

// #[derive(Debug, Deserialize)]
struct SeqRecord {
    seq: Option<String>,
    annotations: Annotations,
    id: Option<String>,
    name: Option<String>,
    description: Option<String>,
    features: Vec<SeqFeature>,
}

// #[derive(Debug, Deserialize)]
struct Annotations {
    molecule_type: Option<String>,
    topology: Option<String>,
    data_file_division: Option<String>,
    date: Option<NaiveDate>,
}

// #[derive(Debug, Deserialize)]
struct SeqFeature {
    location: SimpleLocation,
    feature_type: String,
    qualifiers: Qualifiers,
}

// #[derive(Debug, Deserialize)]
struct SimpleLocation {
    start: usize,
    end: usize,
    strand: Option<i8>,
}

// #[derive(Debug, Deserialize)]
struct Qualifiers {
    parts: Option<String>,
    label: Option<Vec<String>>,
    name: Option<Vec<String>>,
}

#[derive(Debug)]
#[repr(u8)]
// todo: We are observing other packet types: 0x9, 0x3, 0x11, 0x8, 0xd, 0xe, 0x1c
enum PacketType {
    Dna = 0x00,
    Primers = 0x05,
    Notes = 0x06,
    Features = 0x0A,
    Unknown = 0x99, // Placeholder for encountering one we don't recognize
}

impl PacketType {
    pub fn from_byte(byte: u8) -> io::Result<Self> {
        match byte {
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

fn parse_features_packet(payload: &[u8]) -> io::Result<Vec<Feature>> {
    let xml = str::from_utf8(payload).unwrap();
    // let features: Features = from_str(xml).unwrap();

    // Process features data
    Ok(())
}

// #[derive(Debug, Deserialize)]
struct Notes {
    type_: Option<String>,
    last_modified: Option<String>,
    accession_number: Option<String>,
    comments: Option<String>,
}

// #[derive(Debug, Deserialize)]
struct Features {
    // Define fields based on the expected XML structure
}

pub fn import_snapgene(path: &Path) -> io::Result<(SnapgeneData)> {
    let mut file = File::open(path)?;
    parse(&mut file)
}
