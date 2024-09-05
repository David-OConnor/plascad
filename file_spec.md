# The PlasCAD file format

This document describes the PlasCAD file format: a compact way of storing DNA sequences, features, primers, metadata,
and related information. It's a binary format divided into discrete packets, and uses the `.pcad` file extension.
Code for implementing can be found in [pcad.rs](src/file_io/pcad.rs).

Most data structures use [Bincode](https://docs.rs/bincode/latest/bincode/) library. This is convenient for this program's
purposes, but makes external interoperability more challenging.

Byte order is big endian.


## Components
The starting two bytes of a PlasCAD file are always `0xca`, `0xfe`.

The remaining bytes are divided into adjacent packets. Packets can be found in any order.


## Packet structure
- Byte 0: Always `0x11`
- Bytes 1-4: A 32-bit unsigned integer of payload size, in bytes.
- Byte 5: An 8-bit unsigned integer corresponding to the packet's type.
- Bytes 6-end: The payload; how this is encoded depends on packet type.

## Packet types

### Sequence: 0
Contains a DNA sequence.

Bytes 0-3: A 32-bit unsigned integer of sequence length, in nucleotides.
Remaining data: Every two bits is a nucleotide. This packet is always a whole number of bytes; bits in the final byte that
would extend past the sequence length are ignored. Bit assignments for nucleotides is as follows:

- A: `0b00`
- T: `0b01`
- G: `0b10`
- C: `0b11`

#### An example
The sequence `CTGATTTCTG`. This would serialize as follows, using 7 total bytes:
Bytes 0-3: [0, 0, 0, 10], to indicate the sequence length of 10 nucleodies. 

Three additional bytes to encode the sequence; each byte can fit 4 nucleotides:
- Byte 4: `CTGA` `0b11_01_10_00`
- Byte 5: `TTTC` `0b01_01_01_11`
- Byte 6: `TG` `0b01_10_00_00`. 

On the final byte, note the 0-fill on the right; we know not to encode it as `A` due to the
sequence length.

### Features: 1
A bincode serialization of a `Vec<Feature>`

### Primers: 2
A bincode serialization of a `Vec<Primer>`

### Metadata: 3
A bincode serialization of a `Metadata`

### IonConcentrations: 6
A bincode serialization of a `IonConcentrations`

### Portions: 7
A bincode serialization of a `Portions`

### PathLoaded: 10
A bincode serialization of a `Option<PathBuf>`

### Topology: 11
A bincode serialization of a `Topology`



