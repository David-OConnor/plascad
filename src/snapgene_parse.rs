//! Parse SnapGene DNA files.

// This code is translated from [biopython: SnapGeneIO](https://github.com/biopython/biopython/blob/master/Bio/SeqIO/SnapGeneIO.py)

use std::{
    fs::File,
    io::{self, Read, Seek, SeekFrom},
    str,
};

use chrono::NaiveDate;
use regex::Regex;
use serde::Deserialize;
// use serde_xml_rs::from_str;

#[derive(Debug, Deserialize)]
struct SeqRecord {
    seq: Option<String>,
    annotations: Annotations,
    id: Option<String>,
    name: Option<String>,
    description: Option<String>,
    features: Vec<SeqFeature>,
}

#[derive(Debug, Deserialize)]
struct Annotations {
    molecule_type: Option<String>,
    topology: Option<String>,
    data_file_division: Option<String>,
    date: Option<NaiveDate>,
}

#[derive(Debug, Deserialize)]
struct SeqFeature {
    location: SimpleLocation,
    feature_type: String,
    qualifiers: Qualifiers,
}

#[derive(Debug, Deserialize)]
struct SimpleLocation {
    start: usize,
    end: usize,
    strand: Option<i8>,
}

#[derive(Debug, Deserialize)]
struct Qualifiers {
    parts: Option<String>,
    label: Option<Vec<String>>,
    name: Option<Vec<String>>,
}

fn read_bytes<R: Read>(reader: &mut R, length: usize) -> io::Result<Vec<u8>> {
    let mut buffer = vec![0; length];
    reader.read_exact(&mut buffer)?;
    Ok(buffer)
}

fn iterate<R: Read + Seek>(handle: &mut R) -> io::Result<()> {
    loop {
        let packet_type = match read_bytes(handle, 1)? {
            ref v if v.is_empty() => return Ok(()),
            v => v[0],
        };

        let length = read_bytes(handle, 4)?;
        let length = u32::from_be_bytes(length.try_into().unwrap()) as usize;

        let data = read_bytes(handle, length)?;
        parse_packet(packet_type, length, &data)?;
    }
}

fn parse_packet(packet_type: u8, length: usize, data: &[u8]) -> io::Result<()> {
    match packet_type {
        0x00 => parse_dna_packet(length, data),
        0x05 => parse_primers_packet(length, data),
        0x06 => parse_notes_packet(length, data),
        0x0A => parse_features_packet(length, data),
        _ => Ok(()),
    }
}

fn parse_dna_packet(length: usize, data: &[u8]) -> io::Result<()> {
    // Implement parsing logic for DNA packet
    Ok(())
}

fn parse_primers_packet(length: usize, data: &[u8]) -> io::Result<()> {
    // Implement parsing logic for primers packet
    Ok(())
}

fn parse_notes_packet(length: usize, data: &[u8]) -> io::Result<()> {
    let xml = str::from_utf8(data).unwrap();
    let notes: Notes = from_str(xml).unwrap();

    // Process notes data
    Ok(())
}

fn parse_features_packet(length: usize, data: &[u8]) -> io::Result<()> {
    let xml = str::from_utf8(data).unwrap();
    let features: Features = from_str(xml).unwrap();

    // Process features data
    Ok(())
}

#[derive(Debug, Deserialize)]
struct Notes {
    Type: Option<String>,
    LastModified: Option<String>,
    AccessionNumber: Option<String>,
    Comments: Option<String>,
}

#[derive(Debug, Deserialize)]
struct Features {
    // Define fields based on the expected XML structure
}

fn test() -> io::Result<()> {
    let mut file = File::open("path/to/snapgene/file")?;
    iterate(&mut file)
}
