# PlasCAD

[![Crate](https://img.shields.io/crates/v/plascad.svg)](https://crates.io/crates/plascad)
[![Docs](https://docs.rs/plascad/badge.svg)](https://docs.rs/plascad)

Design software for plasmid and primer creation and validation.


## Installation

### Windows and Linux
[Download, unzip, and run](https://github.com/David-OConnor/plascad/releases). 


### Mac
Compile from source by [downloading and installing Rust](https://www.rust-lang.org/tools/install), then running `cargo install plascad` from a CLI.


## Current functionality

### Primer QC and tuning
Evaluates primers for quality in several metrics: melting temperature, GC content, 3' end stability, repeats, and possibility
of forming self-end dimers.

This allows primer length to be automatically or manually changed to optimize these parameters. This is done by marking
one or more primer end as not having a fixed starting point, and providing more-than-expected nucleotides of the matching 
sequence on this end. The end point can then be adjusted to optimize primer qualities.


### Primer generation for SLIC and FastCloning
Given the sequences for an insert, a vector, and insertion point, it will generate primers suitable for SLIC and FastCloning.
It generates primers to amplify the entire vector, and insert primers that contain suitable overlap regions with the vector.


### Sequence viewer
This shows the sequence of interest (as generated from cloning, or manually input) with primers overlayed based on their
match location. It also displays cut sites for common restriction enzymes, and features loaded from a file, or set by the user.


### Circular map
A circular sequence map of the plasmid being edited, with features and other data displayed


### PCR parameter generation
Generates PCR parameters (eg temperature and time for the various steps), based on product length, primer
melting temperature, polymerase type, and other parameters.


### Interop with SnapGene and FASTA
Can read and write FASTA and SnapGene .dna files.


## Why another plasmid editor
We believe that the more tools available for scientists, the better. In particular, my goal is to make
a fast, lightweight program that's as easy to use as possible, without sacrificing functionality. We also added
some functionality we didn't see elsewhere, like automatic primer tuning, and primer generation for SLIC and FastCloning.


## Near-term plans
- A better save and load system
- QCing plasmids for toxic proteins, and various forms of error
- QC primers for problems in context of plasmids. (Eg multiple binding sites)
- Identifying secondary structures, hairpins etc
- Better sequence view and edit functionality, expanding the sequence view.


## Calculations used 
Our primer melting temperature method used is based on [SantaLucia & Hicks (2004)](https://pubmed.ncbi.nlm.nih.gov/15139820/) It uses each pair of adjacent nucleotides in the
primer sequence to estimate entropy and enthalpy values, used in the following calcluation, where $ΔH$ is enthalphy, $ΔS$ is entropy, and $C_T$ is 
primer Molar concentration:

$$ (1000 * ΔH) / (ΔS + R \times ln(\frac{C_T}{4})) - 273.15 $$

The calculation also includes salt correction, derived from BioPython, using concentrations of $K^+$, $Na^+$, $Mg^{2+}$, and dntp concentration. These are provided by the user, or initiated with defaults.

We calculate the following categories of nucleotide repeats:
- Single nucleotides repeated 4 or more times in a row
- Nucleotide pairs repeated 4 or more times in a row
- Any sequence of 3 or more nucleotides that occurs more than once in a sequence

*Primer quality* is a fuzzy metric, and is calculated as a weighted average from other metrics. It's the quantity our tuning algorithm optimizes when adjusting primer length.