//! For reading AB1 trace files. (Applied Biosystem's sequencing)
//! [BioPython docs](https://biopython.org/wiki/ABI_traces)
//!
//! Adapted directly from this [BioPython code](https://github.com/biopython/biopython/blob/master/Bio/SeqIO/AbiIO.py)
//!
//! We are unable to find the official format spec for AB1 files.

use std::{
    collections::HashMap,
    error::Error,
    fs::File,
    io::{self, ErrorKind, Read, Seek, SeekFrom},
    path::Path,
};
use bio::seq::Seq;
use chrono::{NaiveDate, NaiveTime};
use byteorder::{BigEndian, ReadBytesExt};
use std::io::Cursor;
use bio::io::fastq;
use std::cmp;

use na_seq::Seq;

const HEADER_SIZE: usize = 26; // bytes(?)
const DIR_SIZE: usize = 28; // bytes(?)

#[derive(Debug)]
struct Header {
    a: u16,
    b: [u8; 4], // string of 4 bytes
    c: u32,
    d: [u16; 2],
    e: [u32; 3],
}

impl Header {
    pub fn from_bytes(bytes: [u8; HEADER_SIZE]) -> Self {
        Self {
            a: u16::from_be_bytes(bytes[0..2].try_into().unwrap()),
            b: [bytes[2], bytes[3], bytes[4], bytes[5]],
            c: u32::from_be_bytes(bytes[6..10].try_into().unwrap()),
            d: [
                u16::from_be_bytes(bytes[10..12].try_into().unwrap()),
                u16::from_be_bytes(bytes[12..14].try_into().unwrap()),
            ],
            e: [
                u32::from_be_bytes(bytes[14..18].try_into().unwrap()),
                u32::from_be_bytes(bytes[18..22].try_into().unwrap()),
                u32::from_be_bytes(bytes[22..26].try_into().unwrap()),
            ]
        }
    }
}

#[derive(Debug)]
struct Dir {
    a: [u8; 4], // string of 4 bytes
    b: u32,
    c: [u16; 2],
    d: [u32; 4],
}

impl Dir {
    pub fn from_bytes(bytes: [u8; DIR_SIZE]) -> Self {
        Self {
            a: [bytes[0], bytes[1], bytes[2], bytes[3]],
            b: u32::from_be_bytes(bytes[4..8].try_into().unwrap()),
            c: [
                u16::from_be_bytes(bytes[8..10].try_into().unwrap()),
                u16::from_be_bytes(bytes[10..12].try_into().unwrap()),
            ],
            d: [
                u32::from_be_bytes(bytes[12..16].try_into().unwrap()),
                u32::from_be_bytes(bytes[16..20].try_into().unwrap()),
                u32::from_be_bytes(bytes[20..24].try_into().unwrap()),
                u32::from_be_bytes(bytes[24..28].try_into().unwrap()),
            ],
        }
    }
}


#[derive(Debug)]
struct SeqRecord {
    id: String,
    name: String,
    description: String,
    sequence: Option<Seq>,
    annotations: HashMap<String, String>,
    phred_quality: Option<Vec<u8>>,
}

#[derive(Debug)]
struct AbiIterator<R: Read + Seek> {
    stream: R,
    trim: bool,
}

impl<R: Read + Seek> AbiIterator<R> {
    pub fn new(mut stream: R, trim: bool) -> Result<Self, Box<dyn Error>> {
        let mut marker = [0u8; 4];
        stream.read_exact(&mut marker)?;
        println!("Marker: {:?}", marker);
        if &marker != b"ABIF" {
            return Err(format!("File should start with ABIF, not {:?}", marker).into());
        }
        Ok(Self { stream, trim })
    }

    pub fn next(&mut self) -> Result<Option<SeqRecord>, Box<dyn Error>> {
        let mut header_data = [0; HEADER_SIZE];

        if self.stream.read(&mut header_data)? == 0 {
            println!("EOF");
            return Ok(None); // End of file
        }

        let header = Header::from_bytes(header_data);

        println!("Header: {:?}", header);

        // Parse header and directories (simplified for brevity)
        // let raw = HashMap::new(); // Placeholder for raw data
        let sequence = None; // Placeholder for sequence data
        let phred_quality = None; // Placeholder for quality values
        let annotations = HashMap::new(); // Placeholder for annotations

        let record = SeqRecord {
            id: "<unknown id>".to_string(),
            name: "example".to_string(),
            description: "example description".to_string(),
            sequence,
            annotations,
            phred_quality,
        };
        Ok(Some(record))
    }
}

// Helper function to parse ABI tags
// fn parse_abi_tag(data: &[u8]) -> Result<(String, String), Box<dyn Error>> {
fn parse_abi_tag(data: &[u8]) -> Result<(String, String), Box<dyn Error>> {
    let tag_name = String::from_utf8_lossy(&data[0..4]).to_string();
    let tag_number = u32::from_be_bytes(data[4..8].try_into()?);
    Ok((tag_name, tag_number.to_string()))
}


fn abi_trim(seq_record: &fastq::Record) -> fastq::Record {
    // Richard Mott's modified trimming algorithm.

    let segment = 20; // Minimum sequence length
    let cutoff = 0.05; // Default cutoff value for calculating base score

    // If the length of the sequence is less than or equal to the segment size, return as is.
    if seq_record.seq().len() <= segment {
        return seq_record.clone();
    }

    // Calculate base scores from quality values.
    let score_list: Vec<f64> = seq_record
        .qual()
        .iter()
        .map(|&qual| cutoff - 10f64.powf((qual as f64) / -10.0))
        .collect();

    // Calculate cumulative score, initialize with zero.
    let mut cumulative_scores: Vec<f64> = vec![0.0];
    let mut trim_start = 0;
    let mut start_flag = false;

    for i in 1..score_list.len() {
        let score = cumulative_scores[i - 1] + score_list[i];
        if score < 0.0 {
            cumulative_scores.push(0.0);
        } else {
            cumulative_scores.push(score);
            if !start_flag {
                // Set trim_start when cumulative score is first greater than zero
                trim_start = i;
                start_flag = true;
            }
        }
    }

    // Find the index of the highest cumulative score to mark the end of the trimming segment
    let trim_finish = cumulative_scores
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(idx, _)| idx)
        .unwrap_or(0);

    // Extract the trimmed sequence
    let trimmed_seq = &seq_record.seq()[trim_start..trim_finish];
    let trimmed_qual = &seq_record.qual()[trim_start..trim_finish];

    // Create a new trimmed record and return it
    fastq::Record::with_attrs(seq_record.id(), None, trimmed_seq, trimmed_qual)
}

fn parse_tag_data(elem_code: u8, elem_num: usize, raw_data: &[u8]) -> Option<String> {
    let byte_fmt: HashMap<u8, &str> = vec![
        (2, "u8"),
        (10, "date"),
        (11, "time"),
        (13, "bool"),
        (18, "slice_after"),
        (19, "slice_before"),
    ]
        .into_iter()
        .collect();

    if let Some(&fmt) = byte_fmt.get(&elem_code) {
        let mut cursor = Cursor::new(raw_data);
        match fmt {
            "u8" => {
                if elem_num == 1 {
                    Some(cursor.read_u8().ok()?.to_string())
                } else {
                    Some(
                        (0..elem_num)
                            .filter_map(|_| cursor.read_u8().ok())
                            .map(|v| v.to_string())
                            .collect::<Vec<_>>()
                            .join(", "),
                    )
                }
            }
            "date" => {
                let year = cursor.read_u16::<BigEndian>().ok()?;
                let month = cursor.read_u8().ok()?;
                let day = cursor.read_u8().ok()?;
                Some(NaiveDate::from_ymd_opt(year as i32, month as u32, day as u32)?.to_string())
            }
            "time" => {
                let hour = cursor.read_u8().ok()?;
                let minute = cursor.read_u8().ok()?;
                let second = cursor.read_u8().ok()?;
                Some(NaiveTime::from_hms_opt(hour as u32, minute as u32, second as u32)?.to_string())
            }
            "bool" => Some((cursor.read_u8().ok()? != 0).to_string()),
            "slice_after" => Some(
                raw_data[1..]
                    .iter()
                    .map(|&v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(", "),
            ),
            "slice_before" => Some(
                raw_data[..raw_data.len() - 1]
                    .iter()
                    .map(|&v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(", "),
            ),
            _ => None,
        }
    } else {
        None
    }
}

fn abi_parse_header(header: &[u32], stream: &mut File) -> io::Result<Vec<(String, u32, Option<String>)>> {
    let head_elem_size = header[4] as usize;
    let head_elem_num = header[5] as usize;
    let head_offset = header[7] as u64;

    let mut results = Vec::new();

    for index in 0..head_elem_num {
        let start = head_offset + (index * head_elem_size) as u64;
        stream.seek(SeekFrom::Start(start))?;

        let mut buffer = vec![0u8; 28]; // Assuming _DIRFMT is 28 bytes long
        stream.read_exact(&mut buffer)?;

        let mut cursor = Cursor::new(buffer);
        let tag_name = read_string(&mut cursor, 4)?;
        let tag_number = cursor.read_u32::<BigEndian>()?;
        let elem_code = cursor.read_u16::<BigEndian>()? as u8;
        let _elem_size = cursor.read_u16::<BigEndian>()?;
        let elem_num = cursor.read_u32::<BigEndian>()?;
        let data_size = cursor.read_u32::<BigEndian>()?;
        let data_offset = cursor.read_u32::<BigEndian>()?;
        let tag_offset = cursor.read_u32::<BigEndian>()?;

        let data_offset = if data_size <= 4 {
            tag_offset + 20
        } else {
            data_offset
        } as u64;

        stream.seek(SeekFrom::Start(data_offset))?;
        let mut data = vec![0u8; data_size as usize];
        stream.read_exact(&mut data)?;

        let parsed_data = parse_tag_data(elem_code, elem_num as usize, &data);
        results.push((tag_name, tag_number, parsed_data));
    }

    Ok(results)
}

fn read_string<R: Read>(reader: &mut R, length: usize) -> io::Result<String> {
    let mut buffer = vec![0; length];
    reader.read_exact(&mut buffer)?;
    Ok(String::from_utf8_lossy(&buffer).trim_end_matches(char::from(0)).to_string())
}


/// Read a file in the GenBank format.
/// [Rust docs ref of fields](https://docs.rs/gb-io/latest/gb_io/seq/struct.Seq.html)
// pub fn import_ab1(path: &Path) -> io::Result<()> {
pub fn import_ab1(path: &Path) -> Result<(), Box<dyn Error>> {
    let file = File::open(path)?;
    let mut iterator = AbiIterator::new(file, false)?;
    while let Some(record) = iterator.next()? {
        println!("{:?}", record);
    }
    Ok(())
}
