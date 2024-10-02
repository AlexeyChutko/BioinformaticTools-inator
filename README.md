## **BioinformaticTools-inator**
 BioinformaticTools-inator is a toolkit for working with DNA and RNA sequences and filtering fastq sequences. This toolkit includes functions such as `run_dna_rna_tools` and `filter_fastq`.
 
### Installation
- To run BioinformaticTools-inator you need to have Python installed
- Get the code:
```
git clone git@github.com:AlexeyChutko/BioinformaticTools-inator.git
git branch -all
git checkout hw4
```

### Usage
1. run_dna_rna_tools
This toolkit works with DNA and RNA sequences by performing operations such as transcription, complementarity, reverse and nucleotide_frequency calculation
Example:
```
from Tools import run_dna_rna_tools
sequences = 'ATG'
procedure = 'transcribe
print(run_dna_rna_tools(sequences, procedure))
```
2. filter_fastq
This utile filters fastq sequences witht heir length, GC composition and quality (based on phred33)
Example:
```
from Tools import filter_fastq
seqs = '@SRX079873': ('ACAGCA', 'FGGGFG')
gc_bounds = (20, 80)
length_bounds = (30, 70)
quality_threshold = 90
print(filter_fastq(seqs, gc_bounds, length_bounds, quality_threshold))
```
