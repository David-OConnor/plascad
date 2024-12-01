//! For reading AB1 trace files. (Applied Biosystem's sequencing)
//! [BioPython docs](https://biopython.org/wiki/ABI_traces)
//! [BioPython code](https://github.com/biopython/biopython/blob/master/Bio/SeqIO/AbiIO.py)

use std::{
    collections::HashMap,
    error::Error,
    fs::File,
    io::{self, ErrorKind, Read, Seek, SeekFrom},
    path::Path,
};

use na_seq::Seq;
// use byteorder::{BigEndian, ReadBytesExt};

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
        if &marker != b"ABIF" {
            return Err(format!("File should start with ABIF, not {:?}", marker).into());
        }
        Ok(Self { stream, trim })
    }

    pub fn next(&mut self) -> Result<Option<SeqRecord>, Box<dyn Error>> {
        let mut header = [0u8; 24];
        if self.stream.read(&mut header)? == 0 {
            return Ok(None); // End of file
        }
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

/// Read a file in the GenBank format.
/// [Rust docs ref of fields](https://docs.rs/gb-io/latest/gb_io/seq/struct.Seq.html)
// pub fn import_ab1(path: &Path) -> io::Result<()> {
pub fn import_ab1(path: &Path) -> Result<(), Box<dyn Error>> {
    let file = File::open("example.ab1")?;
    let mut iterator = AbiIterator::new(file, false)?;
    while let Some(record) = iterator.next()? {
        println!("{:?}", record);
    }
    Ok(())
}
