# Plasmid tools

This program is a set of tools for assisting with plasmid creation and assessment.

It is a work in progress, and the current feature set is limited. It is currently similar to AmplifX's 
primer validation tool. Expect a gradual increase in new features.


## Current functionality

### Primer QC and tuning
Evaluates primers for quality in several metrics: melting temperature, GC content, 3' end stability, complexity, and
possibility of forming self-end dimers.

Allows primer length to be automatically or manually tuned to optimize these parameters.

(todo: Add details on the calculations and assumptions used)


### Primer generation for SLIC and FastCloning
Given the sequences for an insert, a vector, and insertion point, it will generate primers suitable for SLIC and FastCloning.
It generates primers to amplify the entire vector, and insert primers that contain suitable overlap regions with the vector.


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