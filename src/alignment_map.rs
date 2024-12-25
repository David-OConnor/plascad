//! Data structures and related code for BAM and SAM alignment map data.
//! [SAM/BAM spec document, CAO Nov 2024](https://samtools.github.io/hts-specs/SAMv1.pdf)
//!
//! BAM format is little endian.

use std::fs::File;
use std::io;
use std::io::{ErrorKind, Read};
use std::path::Path;
use flate2::read::GzDecoder;

const MAGIC: [u8; 4] = [b'B', b'A', b'M', 1];

#[derive(Debug)]
struct RefSeq {
    pub name: String,
    pub l_ref: u32
}

impl RefSeq {
    /// Deserialize from BAM.
    pub fn from_buf(buf: &[u8]) -> io::Result<Self> {
        let l_name = u32::from_le_bytes(buf[0..4].try_into().unwrap());
        let name_end = 4 + l_name as usize;

        Ok(
            Self {
                // todo: Don't unwrap this string parsing; map the error
                // name_end - 1: Don't parse the trailing NUL char. (In the spec, this is a null-terminated string)
                name: String::from_utf8(buf[4..name_end - 1].to_vec()).unwrap(),
                l_ref: u32::from_le_bytes(buf[name_end..name_end + 4].try_into().unwrap()),
            }
        )
    }

    /// Serialize to BAM.
    pub fn to_buf(&self) -> Vec<u8> {
        Vec::new()
    }
}

#[derive(Debug)]
struct AuxData {
    pub tag: [u8; 2],
    pub val_type: u8,
    /// corresponds to val_type.
    pub value: Vec<u8>, // todo: Placeholder
}

#[derive(Debug)]
struct Alignment {
    pub block_size: u32,
    pub ref_id: i32,
    pub pos: i32,
    pub l_read_name: u8,
    pub mapq: u8,
    pub bin: u16,
    pub n_cigar_op: u16,
    pub flag: u16,
    pub l_seq: u32,
    pub next_ref_id: i32,
    pub next_pos: i32,
    pub tlen: i32,
    pub read_name: String,
    pub cigar: Vec<u32>,
    pub seq: Vec<u8>, // 4-bit encoded
    pub qual: String, // todo: Vec<u8>?
    pub aux_data: Vec<AuxData>,
}

impl Alignment {
    /// Deserialize from BAM.
    pub fn from_buf(buf: &[u8]) -> io::Result<Self> {
        if buf.len() < 36 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Buffer is too small to contain a valid Alignment record",
            ));
        }

        let block_size = u32::from_le_bytes(buf[0..4].try_into().unwrap());
        let ref_id = i32::from_le_bytes(buf[4..8].try_into().unwrap());
        let pos = i32::from_le_bytes(buf[8..12].try_into().unwrap());
        let l_read_name = buf[12];
        let mapq = buf[13];
        let bin = u16::from_le_bytes(buf[14..16].try_into().unwrap());
        let n_cigar_op = u16::from_le_bytes(buf[16..18].try_into().unwrap());
        let flag = u16::from_le_bytes(buf[18..20].try_into().unwrap());
        let l_seq = u32::from_le_bytes(buf[20..24].try_into().unwrap());
        let next_ref_id = i32::from_le_bytes(buf[24..28].try_into().unwrap());
        let next_pos = i32::from_le_bytes(buf[28..32].try_into().unwrap());
        let tlen = i32::from_le_bytes(buf[32..36].try_into().unwrap());

        let mut offset = 36;
        let read_name = String::from_utf8(buf[offset..offset + l_read_name as usize - 1].to_vec())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        offset += l_read_name as usize;

        let cigar_len = n_cigar_op as usize;
        let mut cigar = Vec::new();
        for _ in 0..cigar_len {
            let cigar_op = u32::from_le_bytes(buf[offset..offset + 4].try_into().unwrap());
            cigar.push(cigar_op);
            offset += 4;
        }

        let seq_len = ((l_seq + 1) / 2) as usize; // Each nucleotide is 4 bits
        let seq = buf[offset..offset + seq_len].to_vec();
        offset += seq_len;

        let qual = buf[offset..offset + l_seq as usize].to_vec();
        offset += l_seq as usize;

        let aux_data = Vec::new(); // TODO: Parse auxiliary data from the remaining bytes if applicable

        Ok(Self {
            block_size,
            ref_id,
            pos,
            l_read_name,
            mapq,
            bin,
            n_cigar_op,
            flag,
            l_seq,
            next_ref_id,
            next_pos,
            tlen,
            read_name,
            cigar,
            seq,
            qual: String::new(), // Placeholder until auxiliary data parsing
            aux_data,
        })
    }

    /// Serialize to BAM.
    pub fn to_buf(&self) -> Vec<u8> {
        Vec::new()
    }
}


#[derive(Debug)]
/// Represents BAM or SAM data. See the spec document, section 4.2. This maps directly to the BAM format.
/// Fields, and their types, here and in sub-structs are taken directly from this table.
pub struct AlignmentMap {
    pub l_text: u32,
    pub text: String,
    pub n_ref: u32,
    pub refs: Vec<RefSeq>,
    pub alignments: Vec<Alignment>
}

impl AlignmentMap {
    /// Deserialize from BAM.
    pub fn from_buf(buf: &[u8]) -> io::Result<Self> {
        println!("Buf len: {:?}", buf.len());
        let l_text = u32::from_le_bytes(buf[4..8].try_into().unwrap());
        let text_end = 8 + l_text as usize;
        let n_ref = u32::from_le_bytes(buf[text_end..text_end + 4].try_into().unwrap());
        let refs_start = text_end + 4;

        let mut i = refs_start;
        let mut refs = Vec::new();
        for _ in 0..n_ref {
            let ref_seq = RefSeq::from_buf(buf[i..].try_into().unwrap())?;
            i += 9 + ref_seq.name.len();

            refs.push(ref_seq);
        }

        let mut alignments = Vec::new();

        while i + 1 < buf.len() {
            let alignment = Alignment::from_buf(buf[i..].try_into().unwrap())?;
            i += 4 + alignment.block_size as usize;

            alignments.push(alignment);
        }

        if buf[0..4] != MAGIC {
            return Err(io::Error::new(ErrorKind::InvalidData, "Incorrect BAM magic.".to_string()));
        }

        Ok(
            Self {
                l_text,
                // todo: Map the error; don't unwrap.
                text: String::from_utf8(buf[8..text_end].to_vec()).unwrap(),
                n_ref,
                refs,
                alignments,
            }
        )
    }

    /// Serialize to BAM.
    pub fn to_buf(&self) -> Vec<u8> {
        Vec::new()
    }
}

/// todo: Move to `file_io` A/R.
pub fn import(path: &Path) -> io::Result<AlignmentMap> {
    let file = File::open(path)?;
    let mut decoder = GzDecoder::new(file);

    println!("Decoding file...");

    let mut decompressed_data = Vec::new();
    decoder.read_to_end(&mut decompressed_data)?;

    println!("Decode complete. Reading...");

    AlignmentMap::from_buf(&decompressed_data)
}