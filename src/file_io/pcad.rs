//! This module contains code for saving and loading in our own PCAD format.
//!
//! This is a binary format that uses packets for each of several message types. This sytem
//! should have better backwards compatibility than raw serialization and deserialization using
//! Bincode or similar.
//!
//! Byte encoding is big endian.

use std::{io, io::ErrorKind};

use bincode::config;
use num_enum::TryFromPrimitive;

use crate::file_io::{
    save::{deser_seq_bin, serialize_seq_bin, StateToSave},
    GenericData,
};

const START_BYTES: [u8; 2] = [0x1f, 0xb2]; // Arbitrary, used as a sanity check.
const PACKET_START: u8 = 0x11;
const PACKET_OVERHEAD: usize = 6; // packet start, packet type, message size.

#[repr(u8)]
#[derive(Clone, Copy, PartialEq, TryFromPrimitive)]
pub enum PacketType {
    Sequence = 0,
    Features = 1,
    Primers = 2,
    Metadata = 3,
    IonConcentrations = 6,
    Portions = 7,
    PathLoaded = 10,
    Topology = 11,
}

/// Byte 0: Standard packet start. Bytes 1-4: u32 of payload len. Bytes 5[..]: Payload.
pub struct Packet {
    type_: PacketType,
    payload: Vec<u8>,
}

impl Packet {
    /// Note: Bytes includes this payload, and potentially until the end of the entire file data.
    /// We use the length bytes to know when to stop reading.
    pub fn from_bytes(bytes: &[u8]) -> io::Result<Self> {
        if bytes.len() < PACKET_OVERHEAD {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Packet must be at least 2 bytes",
            ));
        }

        if bytes[0] != PACKET_START {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Invalid packet start byte in PCAD file.",
            ));
        }

        // todo: This may not be necessary, due to the vec being passed.
        let payload_size = u32::from_be_bytes(bytes[1..5].try_into().unwrap()) as usize;
        if bytes.len() < PACKET_OVERHEAD + payload_size {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Remaining payload is too short based on read packet len in PCAD file.",
            ));
        }

        Ok(Self {
            type_: bytes[5].try_into().map_err(|_| {
                io::Error::new(ErrorKind::InvalidData, "Invalid packet type received")
            })?,
            payload: bytes[PACKET_OVERHEAD..PACKET_OVERHEAD + payload_size].to_vec(), // todo: This essentially clones? Not ideal.
        })
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = vec![PACKET_START];

        let len = self.payload.len() as u32;
        result.extend(&len.to_be_bytes());

        result.push(self.type_ as u8);
        result[6..].copy_from_slice(&self.payload); // todo: Can we avoid this clone?

        result
    }
}

impl StateToSave {
    // Move: Add our packet-related format code to a separate file A/R. eg pcad.rs
    /// Serialize as bytes, eg for file loading. We  use a packet-based format.
    pub fn to_bytes(&self) -> Vec<u8> {
        let cfg = config::standard();
        let mut result = Vec::new();

        result.extend(&START_BYTES);

        // Note: The order we add these packets in doesn't make a difference for loading.

        let seq_packet = Packet {
            type_: PacketType::Sequence,
            payload: serialize_seq_bin(&self.generic.seq),
        };

        let features_packet = Packet {
            type_: PacketType::Features,
            payload: bincode::encode_to_vec(&self.generic.features, cfg).unwrap(),
        };

        let primers_packet = Packet {
            type_: PacketType::Primers,
            payload: bincode::encode_to_vec(&self.generic.primers, cfg).unwrap(),
        };

        let metadata_packet = Packet {
            type_: PacketType::Metadata,
            payload: bincode::encode_to_vec(&self.generic.metadata, cfg).unwrap(),
        };

        let topology_packet = Packet {
            type_: PacketType::Topology,
            payload: bincode::encode_to_vec(&self.generic.topology, cfg).unwrap(),
        };

        let ion_concentrations_packet = Packet {
            type_: PacketType::IonConcentrations,
            payload: bincode::encode_to_vec(&self.ion_concentrations, cfg).unwrap(),
        };

        let portions_packet = Packet {
            type_: PacketType::Portions,
            payload: bincode::encode_to_vec(&self.portions, cfg).unwrap(),
        };

        let path_loaded_packet = Packet {
            type_: PacketType::PathLoaded,
            payload: bincode::encode_to_vec(&self.path_loaded, cfg).unwrap(),
        };

        result.extend(&seq_packet.to_bytes());
        result.extend(&features_packet.to_bytes());
        result.extend(&primers_packet.to_bytes());
        result.extend(&metadata_packet.to_bytes());
        result.extend(&topology_packet.to_bytes());

        result.extend(&ion_concentrations_packet.to_bytes());
        result.extend(&portions_packet.to_bytes());
        result.extend(&path_loaded_packet.to_bytes());

        result
    }

    pub fn from_bytes(bytes: &[u8]) -> io::Result<Self> {
        if bytes[0..2] != START_BYTES {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Invalid start bytes in PCAD file.",
            ));
        }

        const DECODE_ERR_MSG: &str = "Error parsing a packet component";

        let cfg = config::standard();
        let mut result = StateToSave::default();

        // for(i, byte) in bytes[2..].iter().enumerate() {
        let mut i = 2;
        loop {
            if i > bytes.len() {
                break; // End of the packet.
            }

            let bytes_remaining = &bytes[i..];
            let packet = Packet::from_bytes(&bytes_remaining)?;
            i += PACKET_OVERHEAD + packet.payload.len();

            // Now, add packet data to our result A/R.
            match packet.type_ {
                PacketType::Sequence => {
                    result.generic.seq = deser_seq_bin(&packet.payload)
                        .map_err(|_| io::Error::new(ErrorKind::InvalidData, DECODE_ERR_MSG))?;
                }
                PacketType::Features => {
                    result.generic.features = bincode::decode_from_slice(&packet.payload, cfg)
                        .map_err(|_| io::Error::new(ErrorKind::InvalidData, DECODE_ERR_MSG))?.0;
                }
                PacketType::Primers => {
                    result.generic.primers = bincode::decode_from_slice(&packet.payload, cfg)
                        .map_err(|_| io::Error::new(ErrorKind::InvalidData, DECODE_ERR_MSG))?.0;
                }
                PacketType::Metadata => {
                    result.generic.metadata = bincode::decode_from_slice(&packet.payload, cfg)
                        .map_err(|_| io::Error::new(ErrorKind::InvalidData, DECODE_ERR_MSG))?.0;
                }
                PacketType::Topology => {
                    result.generic.topology = bincode::decode_from_slice(&packet.payload, cfg)
                        .map_err(|_| io::Error::new(ErrorKind::InvalidData, DECODE_ERR_MSG))?.0;
                }
                PacketType::IonConcentrations => {
                    result.ion_concentrations = bincode::decode_from_slice(&packet.payload, cfg)
                        .map_err(|_| io::Error::new(ErrorKind::InvalidData, DECODE_ERR_MSG))?.0;
                }
                PacketType::Portions => {
                    result.portions = bincode::decode_from_slice(&packet.payload, cfg)
                        .map_err(|_| io::Error::new(ErrorKind::InvalidData, DECODE_ERR_MSG))?.0;
                }
                PacketType::PathLoaded => {
                    result.path_loaded = bincode::decode_from_slice(&packet.payload, cfg)
                        .map_err(|_| io::Error::new(ErrorKind::InvalidData, DECODE_ERR_MSG))?.0;
                }
            }
        }

        Ok(result)
    }
}
