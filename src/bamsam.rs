//! Data structures and related code for BAM and SAM alignment map data.
//! [SAM/BAM spec document, CAO Nov 2024](https://samtools.github.io/hts-specs/SAMv1.pdf)
//!
//! BAM format is little endian.

use std::io;

struct RefSeq {
    pub l_name: u32,
    pub name: String,
    pub l_ref: u32
}

impl RefSeq {
    /// Deserialize from BAM.
    pub fn from_buf(buf: &[u8]) -> io::Result<Self> {
        let l_name = u32::from_le_bytes(buf[0..4]);
        let name_end = 4 + l_name as usize;

        Ok(
            Self {
                l_name,
                name: String::from_utf8(buf[4..name_end]),
                l_ref: u32::from_le_bytes(buf[name_end..name_end + 4]),
            }
        )
    }

    /// Serialize to BAM.
    pub fn to_buf(&self) -> Vec<u8> {
        Vec::new();
    }
}

struct AuxData {
    pub tag: [u8; 2],
    pub val_type: u8,
    /// corresponds to val_type.
    pub value: Vec<u8>, // todo: Placeholder
}

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
    pub read_name: i32,
    pub cigar: String, // todo: Vec<u8>?
    pub seq: Vec<u8>, // 4-bit encoded
    pub qual: String, // todo: Vec<u8>?
    pub aux_data: Vec<AuxData>,
}

impl Alignment {
    /// Deserialize from BAM.
    pub fn from_buf(buf: &[u8]) -> io::Result<Self> {
        Ok(
            Self {

            }
        )
    }

    /// Serialize to BAM.
    pub fn to_buf(&self) -> Vec<u8> {
        Vec::new();
    }
}


/// See the spec document, section 4.2. This maps directly to the BAM format.
/// Fields, and their types, here and in sub-structs are taken directly from this table.
struct SamBamData {
    pub magic: [u8; 4],
    pub l_text: u32,
    pub text: String,
    pub n_ref: u32,
    pub refs: Vec<RefSeq>,
    pub alignments: Vec<Alignment>
}

impl SamBamData {
    /// Deserialize from BAM.
    pub fn from_buf(buf: &[u8]) -> io::Result<Self> {
        let l_text = u32::from_le_bytes(buf[4..8]);
        let text_end = 8 + l_text as usize;
        let refs_start = text_end + 4;
        let alignments_start = 0; // todo

        let mut i = refs_start;
        let mut refs = Vec::new();
        while i + 1 < buf.len() {
            let ref_seq = RefSeq::from_buf(buf[refs_start..]);
            refs.push(ref_seq);
            i += 8 + ref_seq.l_name as usize; // todo: What is the l_ref field for? Is that a more direct way?
        }

        let mut alignments = Vec::new();

        Ok(
            Self {
                magic: buf[0..4],
                l_text,
                text: String::from_utf8(buf[8..text_end].to_vec()),
                n_ref: u32::from_le_bytes(buf[text_end..text_end + 4]),
                refs,
                alignments,
            }
        )
    }

    /// Serialize to BAM.
    pub fn to_buf(&self) -> Vec<u8> {
        Vec::new();
    }
}
