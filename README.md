# Plasmid tools

This program is a set of tools for assisting with plasmid creation and assessment.

It is a work in progress, and the current feature set is limited. It is currently similar to AmplifX's 
primer validation tool. Expect a gradual increase in new features.

[Download for Windows, and Linux](https://github.com/David-OConnor/plasmid_tools/releases). (Unzip, and run). If using Mac, you can compile from source by [downloading Rust](https://www.rust-lang.org/tools/install), and running `cargo install plasmid_tools` in the CLI.

## Current functionality

### Primer QC and tuning
Evaluates primers for quality in several metrics: melting temperature, GC content, 3' end stability, repeats, and possibility of forming self-end dimers.

This allows primer length to be automatically or manually changed to optimize these parameters. This is done by marking one or more primer end as not having a fixed starting point, and providing more-than-expected nucleotides of the matching sequence on this end. The end point can then be adjusted to optimize primer qualities.


### Primer generation for SLIC and FastCloning
Given the sequences for an insert, a vector, and insertion point, it will generate primers suitable for SLIC and FastCloning.
It generates primers to amplify the entire vector, and insert primers that contain suitable overlap regions with the vector.


### A simple sequence viewer
This shows the sequence of interest (as generated from cloning, or manually input) with primers overlayed based on their match location.


### PCR parameter generation
Generates PCR parameters (eg temperature and time for the various steps), based on product length, primer
melting temperature, polymerase type, and other parameters.

(todo: Add details on the calculations and assumptions used)


## Near-term plans
- Visualization of primers along sequences
- A better save and load system
- Reading and writing FASTA, FASTQ, and SnapGene (.dna) files
- QCing plasmids for toxic proteins, and various forms of error
- QC primers for problems in context of plasmids. (Eg multiple binding sites)
- Plasmid view and edit functionality, expanding the sequence view.


## Calculations used 
Our primer melting temperature method used is based on [BioPython's](https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html). This is, in turn, derived from 
[SantaLucia & Hicks (2004)](https://pubmed.ncbi.nlm.nih.gov/15139820/) It uses the neighboring nucleotides of each nucleotide in the sequence to estimate 
entropy and enthalpy values, used in the following calcluation, where ΔH is enthalphy, ΔS is entropy, and $C_T$ is 
primer Molar concentration:

$$ (1000 * ΔH) / (ΔS + R \times ln(\frac{C_T}{4})) - 273.15 $$

It includes salt correction, also derived from BioPython, using concentrations of K+, Na+, Mg2+, and dntp concentration, provided by the user, or initiated with defaults.

We calculate the following sorts of nucleotide repeats:
- Single nucleotides repeated 4 or more times in a row
- Nucleotide pairs repeated 4 or more times in a row
- Any sequence of 3 nucleotides or more that occurs more than once in a sequence

Primer quality is a fuzzy metric, and is calculated as a weighted average from other metrics. It's the quantity our tuning algorithm optimizes when adjusting primer length.