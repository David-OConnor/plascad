//! Data structures and related code for BAM and SAM alignment map data.
//! [SAM/BAM spec document, CAO Nov 2024](https://samtools.github.io/hts-specs/SAMv1.pdf)
//!
//! BAM format is little endian.

use std::{
    fs::File,
    io,
    io::{BufRead, ErrorKind, Read},
    path::Path,
};

use bgzip::{index::BGZFIndex, read::IndexedBGZFReader, BGZFError, BGZFReader};
use flate2::read::{GzDecoder, MultiGzDecoder};

/// Syntax helper for parsing multi-byte fields into primitives.
///
/// Example: `parse_le!(bytes, i32, 5..9);`
#[macro_export]
macro_rules! parse_le {
    ($bytes:expr, $t:ty, $range:expr) => {{
        <$t>::from_le_bytes($bytes[$range].try_into().unwrap())
    }};
}

const MAGIC: [u8; 4] = [b'B', b'A', b'M', 1];

#[derive(Debug)]
struct RefSeq {
    pub name: String,
    pub l_ref: u32,
}

impl RefSeq {
    /// Deserialize from BAM.
    pub fn from_buf(buf: &[u8]) -> io::Result<Self> {
        let l_name = parse_le!(buf, u32, 0..4);
        let name_end = 4 + l_name as usize;

        Ok(Self {
            // todo: Don't unwrap this string parsing; map the error
            // name_end - 1: Don't parse the trailing NUL char. (In the spec, this is a null-terminated string)
            name: String::from_utf8(buf[4..name_end - 1].to_vec()).unwrap(),
            l_ref: parse_le!(buf, u32, name_end..name_end + 4),
        })
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
    /// Size of the entire alignment packet, *except the 4-byte block size*. Includes aux data.
    pub block_size: u32,
    pub ref_id: i32,
    pub pos: i32,
    pub mapq: u8,
    pub bin: u16,
    pub n_cigar_op: u16,
    pub flag: u16,
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

        let block_size = parse_le!(buf, u32, 0..4);

        println!("Block size: {:?}", block_size);

        let ref_id = parse_le!(buf, i32, 4..8);
        let pos = parse_le!(buf, i32, 8..12);

        let l_read_name = buf[12] as usize;
        let mapq = buf[13];

        let bin = parse_le!(buf, u16, 14..16);
        let n_cigar_op = parse_le!(buf, u16, 16..18);
        let flag = parse_le!(buf, u16, 18..20);

        let l_seq = parse_le!(buf, u32, 20..24) as usize;

        let next_ref_id = parse_le!(buf, i32, 24..28);
        let next_pos = parse_le!(buf, i32, 28..32);
        let tlen = parse_le!(buf, i32, 32..36);

        let mut i = 36;
        // -1: Ommit the trailing null.
        println!("Len read name: {:?}", l_read_name);
        let read_name = String::from_utf8(buf[i..i + l_read_name - 1].to_vec())
            .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?;
        i += l_read_name;

        let cigar_len = n_cigar_op as usize;
        let mut cigar = Vec::new();
        for _ in 0..cigar_len {
            let cigar_op = parse_le!(buf, u32, i..i + 4);
            cigar.push(cigar_op);
            i += 4;
        }

        // todo: Figure out how this is actually handled, and set the struct type A/R.
        let seq_len = (l_seq + 1) / 2; // Each nucleotide is 4 bits
        let seq = buf[i..i + seq_len].to_vec();
        i += seq_len;

        let qual = buf[i..i + l_seq].to_vec();
        i += l_seq;

        let aux_data = Vec::new(); // todo: Parse auxiliary data from the remaining bytes.

        Ok(Self {
            block_size,
            ref_id,
            pos,
            mapq,
            bin,
            n_cigar_op,
            flag,
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
    pub fn to_buf_(&self) -> Vec<u8> {
        // todo: If desired.
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
    pub alignments: Vec<Alignment>,
}

impl AlignmentMap {
    /// Deserialize an alignment record from a buffer containing BAM data.
    pub fn from_buf(buf: &[u8]) -> io::Result<Self> {
        if buf.len() < 12 {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Buffer is too small to contain a valid header",
            ));
        }

        let l_text = parse_le!(buf, u32, 4..8);
        let text_end = 8 + l_text as usize;
        let n_ref = parse_le!(buf, u32, text_end..text_end + 4);
        let refs_start = text_end + 4;

        let mut i = refs_start;
        let mut refs = Vec::new();
        for _ in ..n_ref {
            let ref_seq = RefSeq::from_buf(buf[i..].try_into().unwrap())?;
            i += 9 + ref_seq.name.len();

            refs.push(ref_seq);
        }

        let mut alignments = Vec::new();

        while i + 1 < buf.len() {
            let alignment = Alignment::from_buf(buf[i..].try_into().unwrap())?;

            println!("\n\n Alignment: {:?}", alignment);

            i += 4 + alignment.block_size as usize;

            alignments.push(alignment);
        }

        if buf[0..4] != MAGIC {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Incorrect BAM magic.".to_string(),
            ));
        }

        Ok(Self {
            l_text,
            // todo: Map the error; don't unwrap.
            text: String::from_utf8(buf[8..text_end].to_vec()).unwrap(),
            n_ref,
            refs,
            alignments,
        })
    }

    /// Serialize to BAM.
    pub fn to_buf(&self) -> Vec<u8> {
        Vec::new()
    }
}

pub fn import(path: &Path) -> io::Result<AlignmentMap> {
    let file = File::open(path)?;
    let mut decoder = GzDecoder::new(file);
    // let mut decoder = MultiGzDecoder::new(file);

    println!("Decoding file...");

    let mut buf = Vec::new();
    // let mut buf = [0; 80_000];
    decoder.read_to_end(&mut buf)?;

    println!("Decode complete. Reading...");

    let header = AlignmentMap::from_buf(&buf);

    header

    // AlignmentMap::from_buf(&buf)
}

// pub fn import(path: &Path) -> io::Result<AlignmentMap> {
//     let file = File::open(path)?;
//
//     let mut reader =
//         BGZFReader::new(file).unwrap(); // todo: Map error.
//
//     // let index = BGZFIndex::from_reader(reader)?;
//
//     // println!("Index entries: {:?}", index.entries());
//
//     // let mut reader =
//     //     IndexedBGZFReader::new(reader, index).unwrap(); // todo: Map error.
//
//     // println!("Voffset A: {:?}", reader.bgzf_pos());
//
//     reader.bgzf_seek(0).unwrap(); // todo: Map error.
//     // reader.seek(0).unwrap(); // todo: Map error.
//
//     // let mut buf = Vec::new();
//     let mut buf = vec![0; 80_000];
//
//     // reader.read_to_end(&mut buf).unwrap(); // todo: Map error.
//     reader.read(&mut buf).unwrap(); // todo: Map error.
//
//     println!("Voffset: {:?}", reader.bgzf_pos());
//
//     println!("BUf contents: {:?}", &buf[13180..13220]);
//
//     AlignmentMap::from_buf(&buf)
// }

// /// todo: Move to `file_io` A/R.
// /// Note: BAM files are compressed using BGZF, and are therefore split into blocks no larger than 64kb.
// pub fn import(path: &Path) -> io::Result<AlignmentMap> {
//     let mut file = File::open(path)?;
//     let mut decompressed_data = Vec::new();
//
//     // todo: Organize this A/R.
//
//     loop {
//         println!("BAM block...\n");
//         // Read BGZF block header (18 bytes)
//         let mut block_header = [0u8; 18];
//         if let Err(e) = file.read_exact(&mut block_header) {
//             if e.kind() == ErrorKind::UnexpectedEof {
//                 break; // Reached end of file
//             }
//             return Err(e);
//         }
//
//         // Parse block size (last 2 bytes of header, little-endian)
//         let block_size = u16::from_le_bytes([block_header[16], block_header[17]]) as usize;
//
//         println!("Block size: {:?}", block_size);
//
//         if block_size < 18 {
//             return Err(io::Error::new(
//                 ErrorKind::InvalidData,
//                 "Invalid BGZF block size",
//             ));
//         }
//
//         // Read the block data
//         let mut block_data = vec![0u8; block_size - 18];
//         file.read_exact(&mut block_data)?;
//
//         // Decompress the block
//         let mut decoder = GzDecoder::new(&block_data[..]);
//         let mut decompressed_block = Vec::new();
//         decoder.read_to_end(&mut decompressed_block)?;
//
//         // Append to the full decompressed buffer
//         decompressed_data.extend_from_slice(&decompressed_block);
//     }
//
//     AlignmentMap::from_buf(&decompressed_data)
// }

/// Calculate bin given an alignment covering [beg, end) (zero-based, half-closed-half-open)
/// This is adapted from C code in the spec document.
fn reg_to_bin(beg: i32, end: i32) -> i32 {
    let end = end - 1; // Adjust end to be inclusive
    if beg >> 14 == end >> 14 {
        return ((1 << 15) - 1) / 7 + (beg >> 14);
    }
    if beg >> 17 == end >> 17 {
        return ((1 << 12) - 1) / 7 + (beg >> 17);
    }
    if beg >> 20 == end >> 20 {
        return ((1 << 9) - 1) / 7 + (beg >> 20);
    }
    if beg >> 23 == end >> 23 {
        return ((1 << 6) - 1) / 7 + (beg >> 23);
    }
    if beg >> 26 == end >> 26 {
        return ((1 << 3) - 1) / 7 + (beg >> 26);
    }
    0
}

/// Calculate the list of bins that may overlap with region [beg, end) (zero-based)
/// This is adapted from C code in the spec document.
fn reg_to_bins(beg: i32, end: i32, list: &mut Vec<u16>) -> usize {
    let end = end - 1; // Adjust end to be inclusive
    list.push(0);
    for k in (1 + (beg >> 26))..=(1 + (end >> 26)) {
        list.push(k as u16);
    }
    for k in (9 + (beg >> 23))..=(9 + (end >> 23)) {
        list.push(k as u16);
    }
    for k in (73 + (beg >> 20))..=(73 + (end >> 20)) {
        list.push(k as u16);
    }
    for k in (585 + (beg >> 17))..=(585 + (end >> 17)) {
        list.push(k as u16);
    }
    for k in (4681 + (beg >> 14))..=(4681 + (end >> 14)) {
        list.push(k as u16);
    }
    list.len()
}
