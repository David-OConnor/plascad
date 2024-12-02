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
// use bio::seq::Seq;
use chrono::{NaiveDate, NaiveTime};
use byteorder::{BigEndian, ReadBytesExt};
use std::io::Cursor;
use bio::io::fastq;
use std::cmp;

use na_seq::{Nucleotide, Seq, seq_from_str};

const HEADER_SIZE: usize = 26;
const DIR_SIZE: usize = 28;

#[derive(Debug)]
struct Header {
    pub file_version: u16,
    pub tag_name: Seq, // 4 bytes always.
    pub tag_number: u32,
    pub element_type_code: u16,
    pub element_size: u16,
    pub num_elements: usize,
    pub data_size: u32,
    pub data_offset: u32,
}

impl Header {
    pub fn from_bytes(bytes: [u8; HEADER_SIZE]) -> io::Result<Self> {
        let seq_str = std::str::from_utf8(&bytes[2..6]).unwrap().to_owned(); // todo: Handle
        Ok(Self {
            file_version: u16::from_be_bytes(bytes[0..2].try_into().unwrap()),
            tag_name: seq_from_str(&seq_str),
            tag_number: u32::from_be_bytes(bytes[6..10].try_into().unwrap()),
            element_type_code: u16::from_be_bytes(bytes[10..12].try_into().unwrap()),
            element_size: u16::from_be_bytes(bytes[12..14].try_into().unwrap()),
            num_elements: u32::from_be_bytes(bytes[14..18].try_into().unwrap()) as usize,
            data_size: u32::from_be_bytes(bytes[18..22].try_into().unwrap()),
            data_offset: u32::from_be_bytes(bytes[22..26].try_into().unwrap()),
        })
    }
}

#[derive(Debug)]
struct Dir {
    // todo: This breakdown into fields is wrong, but I'm not sure how.
    pub tag_name: String, // 4 bytes
    pub tag_number: u32,
    pub elem_code: u16,
    pub a: u16, // placeholder
    pub num_elements: usize,
    pub data_size: usize,
    pub data_offset: usize,
    pub b: u32,// placeholder
    pub tag_offset: usize,
    // todo: Tag offset??
}

impl Dir {
    pub fn from_bytes(bytes: [u8; DIR_SIZE], tag_offset: usize) -> io::Result<Self> {
        Ok(Self {
            tag_name: std::str::from_utf8(&bytes[..4]).unwrap().to_owned(),// todo: Handle
            tag_number: u32::from_be_bytes(bytes[4..8].try_into().unwrap()),
            elem_code: u16::from_be_bytes(bytes[8..10].try_into().unwrap()),
            a: u16::from_be_bytes(bytes[10..12].try_into().unwrap()),
            num_elements: u32::from_be_bytes(bytes[12..16].try_into().unwrap()) as usize,
            data_size: u32::from_be_bytes(bytes[16..20].try_into().unwrap()) as usize,
            data_offset: u32::from_be_bytes(bytes[20..24].try_into().unwrap()) as usize,
            b: u32::from_be_bytes(bytes[24..28].try_into().unwrap()),
            tag_offset,
        })
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
            return Ok(None); // End of file
        }

        let header = Header::from_bytes(header_data)?;
        println!("Header: {:?}", header);

        for i in 0..header.num_elements as usize {
            // todo: QC data_offset; coming out much too high.
            // Note: Element size should always be DIR_SIZE.
            let start = header.data_offset as usize + i * header.element_size as usize;

            self.stream.seek(SeekFrom::Start(start as u64))?;
            let mut dir_buf = [0; DIR_SIZE];
            if self.stream.read(&mut dir_buf)? == 0 {
                return Ok(None);
            };

            let mut dir = Dir::from_bytes(dir_buf, start)?;

            let key = format!("{}{}", dir.tag_name, dir.tag_number);
            println!("DIR: {:?}, KEY: {:?}", dir, key);

            if dir.data_size <= 4 {
                dir.data_offset = dir.tag_offset + 20;
            }

            self.stream.seek(SeekFrom::Start(dir.data_offset as u64))?;
            let mut tag_buf = vec![0; dir.data_size];
            if self.stream.read(&mut tag_buf)? == 0 {
                return Ok(None)
            };
            // println!("Tag buf: {:?}", tag_buf);

            let tag_data = parse_tag_data(dir.elem_code, dir.num_elements, &tag_buf);

            println!("Tag data: {:?}", tag_data);
        }

        // let raw = HashMap::new(); // Placeholder for raw data

        // for _ in 0..2 {
        //     let key = tag_name + str(tag_number);
        //
        //     // raw[key] = tag_data;
        //
        //     // PBAS2 is base-called sequence, only available in 3530
        //     if key == "PBAS2" {
        //         // seq = tag_data.decode()
        //         // PCON2 is quality values of base-called sequence
        //     } else if key == "PCON2" {
        //         // qual = [ord(val) for val in tag_data.decode()]
        //         // SMPL1 is sample id entered before sequencing run, it must be
        //         //  a string.
        //     } else if key == "SMPL1" {
        //         // sample_id = _get_string_tag(tag_data);
        //     // } else if times.includes(key) {
        //     } else if true {
        //         // times[key] = tag_data;
        //     } else {
        //         if key in _EXTRACT {
        //             annot[_EXTRACT[key]] = tag_data
        //         }
        //     }
        // }

        println!("Header: {:?}", header);

        // Parse header and directories (simplified for brevity)
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

#[derive(Debug)]
enum TagData {
    U8(Vec<u8>),
    U16(Vec<u16>),
}

fn parse_tag_data(elem_code: u16, elem_num: usize, data: &[u8]) -> Option<TagData> {
        //     1: "b",  # byte
        //     2: "s",  # char
        //     3: "H",  # word
        //     4: "h",  # short
        //     5: "i",  # long
        //     6: "2i",  # rational, legacy unsupported
        //     7: "f",  # float
        //     8: "d",  # double
        //     10: "h2B",  # date
        //     11: "4B",  # time
        //     12: "2i2b",  # thumb
        //     13: "B",  # bool
        //     14: "2h",  # point, legacy unsupported
        //     15: "4h",  # rect, legacy unsupported
        //     16: "2i",  # vPoint, legacy unsupported
        //     17: "4i",  # vRect, legacy unsupported
        //     18: "s",  # pString
        //     19: "s",  # cString
        //     20: "2i",  # tag, legacy unsupported

    match elem_code {
        2 => Some(TagData::U8(data.to_vec())),
        4 => {
            let as_u16 = data.chunks_exact(2)
                .map(|chunk| u16::from_le_bytes([chunk[0], chunk[1]])) // Change to `from_be_bytes` for big-endian
                .collect();
            Some(TagData::U16(as_u16))
        },
        _ => {
            // todo: Handle appropriately.
            panic!("Invalid element code in AB1 file: {:?}", elem_code);
        }
    }

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
