## **BioinformaticTools-inator**
 BioinformaticTools-inator is a toolkit for working with DNA and RNA sequences, FASTA files and filtering fastq sequences. This toolkit includes functions such as `run_dna_rna_tools` and `filter_fastq` inside `BioinformaticToolsInator.py` and functions such as `convert_multiline_fasta_to_oneline` and `parse_blast_output` inside `bio_files_processor.py`.
 
### Installation
- To run BioinformaticToolsInator you need to have Python installed
- Get the code:
```
git clone git@github.com:AlexeyChutko/BioinformaticTools-inator.git
git branch -all
git checkout hw15
```

### Usage
#### BioinformaticToolsInator.py
1. **NucleicAcidSequence**  
This toolkit works with DNA and RNA sequences by performing operations such as:  
- `transcribe`: Transcribes a DNA sequence to RNA (replaces thymine 'T' with uracil 'U').
- `complement`: Generates the complement of a DNA or RNA sequence.
- `reverse`:  Reverses the input sequence(s).
- `reverse_complement`: Generates the reverse complement of a DNA or RNA sequence.
- `nucleotide_frequency`: Computes the nucleotide frequency (A, T/U, G, C) for each input sequence.
- `gc_content`: Calculates the GC-content (percentage of G and C nucleotides) of each input sequence.
- `find_start_codons`: Identifies the positions of start codons (AUG for RNA) in each sequence.  

Example for DNASequence:
```
from BioinformaticToolsInator import DNASequence
dna_seq = DNASequence('ATGCATGC')
print(dna_seq.reverse_complement()) # GCATGCAT
```
Example for RNASequence:
```
from BioinformaticToolsInator import RNASequence
rna_seq = RNASequence("AUGCAUGC")
print(rna_seq.find_start_codons()) # [0, 4]
```
2. **AminoAcidSequence**  
This class works with proteinsequences by performing operations such as:  
- `aminoacid_frequency`: Computes the aminoacid frequency for the input sequence.
- `molecular_weight`: Calculates the molecular weight of an amino acid sequence (excluding modifications)

Example for AminoAcidSequence:
```
from BioinformaticToolsInator import AminoAcidSequence
aa = AminoAcidSequence("ACDEFGHIKLMNPQRSTVWY")
print(aa.molecular_weight()) # 2738.1
```
2. **filter_fastq**  
This utile filters fastq sequences with their length, GC composition and quality (based on phred33)  
Filtering Criteria  
- `gc_bounds`: length boundaries within which filtered sequences must be included.
- `length_bounds`: GC-content boundaries (from 0 to 100) within which filtered sequences must be included.
- `quality_threshold`: sequence quality value according toc phred33, below which filtering will not be performed.

Example:
```
from BioinformaticToolsInator import filter_fastq
seqs = {'@SRX079873': ('ACAGCA', 'FGGGFG'), '@SRX079817': ('ATTAGC', 'BFFFFF), '@SRX079858': ('ATGACCCG', 'DCD@@BBC')}
gc_bounds = (20, 80)
length_bounds = (30, 70)
quality_threshold = 90
print(filter_fastq(seqs, gc_bounds, length_bounds, quality_threshold))
```

#### bio_files_processor.py
1. **convert_multiline_fasta_to_oneline**
This utile converts a multi-line FASTA file to a one-line format for each sequence

Example:
```
convert_multiline_fasta_to_oneline("sequences.fasta", "sequences_oneline.fasta")
```
Output:
A new FASTA file where each sequence is written on a single line.

2. **parse_blast_output**
This utile parses BLAST output to extract the best match protein names and saves them to a file

Example:
```
parse_blast_output("blast_results.txt", "parsed_proteins.txt")
```
Output:
A text file containing the sorted list of best match protein names.