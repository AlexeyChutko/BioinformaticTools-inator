from typing import Dict, List, Tuple, Union
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os
import argparse
import logging


class BiologicalSequence(ABC):
    '''
    An abstract base class for representing biological sequences.
    The class provides a general framework for working with biological sequences such as DNA, RNA, and proteins.
    It includes methods for validating the alphabet, obtaining sequence length, and accessing sequence elements.
    '''
    def __init__(self, sequence: str):
        self.sequence = sequence.upper()
        self._validate_alphabet()

    @abstractmethod
    def _validate_alphabet(self):
        pass

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.sequence}')"


class NucleicAcidSequence(BiologicalSequence):
    '''
    Abstract class for nucleic acid sequences (DNA and RNA)
    Can find complement and reversed complement, do reverese,
    calculate nucleotide frequency and GC content of sequences
    '''
    def complement(self):
        '''
        Returns a complementary sequence
        '''
        dna_dict = str.maketrans("ACGTacgt", "TGCAtgca")
        rna_dict = str.maketrans("ACGUacgu", "UGCAugca")
        trans_table = rna_dict if "U" in self.sequence or "u" in self.sequence else dna_dict
        return self.sequence.translate(trans_table)

    def reverse(self):
        '''
        Returns reversed sequence
        '''
        return self.sequence[::-1]

    def reverse_complement(self):
        '''
        Returns reversed complementary sequnce'''
        return self.complement()[::-1]

    def nucleotide_frequency(self):
        '''
        Returns a dictionary with the frequency of each nucleotide in sequence
        '''
        frequency = {}
        for nucl in self.sequence:
            frequency[nucl] = frequency.get(nucl, 0) + 1
        return frequency

    def gc_content(self):
        '''
        Returns a percent of G and C nucleotides ins sequence
        '''
        g_count = self.sequence.count('G')
        c_count = self.sequence.count('C')
        gc_count = g_count + c_count
        return int(round((gc_count / len(self.sequence)) * 100, 0)) if self.sequence else 0

    def _validate_alphabet(self):
        valid_nucleotides = set("ACGTU")
        if not set(self.sequence).issubset(valid_nucleotides):
            raise ValueError("Invalid nucleotide sequence")


class DNASequence(NucleicAcidSequence):
    '''
    Class for DNA sequences.
    Can to transcription of sequences
    '''
    def transcribe(self):
        return self.sequence.replace('T', 'U').replace('t', 'u')


class RNASequence(NucleicAcidSequence):
    '''
    Class for RNA sequences.
    Can find start codons of sequences
    '''
    def find_start_codons(self):
        start_codon = "AUG"
        positions = [i for i in range(len(self.sequence) - 2) if self.sequence[i:i+3] == start_codon]
        return positions


class AminoAcidSequence(BiologicalSequence):
    '''
    Class for amino acid sequences.
    Can calculate molecular weight of sequence and the frequency of each aminoacid in sequence
    '''
    def molecular_weight(self):
        """Calculates the molecular weight of an amino acid sequence (excluding modifications).
        Aminoacid mass was taken from Thermo Fisher Scientific website"""
        amino_acid_weights = {
            'A': 89.1, 'C': 121.2, 'D': 133.1, 'E': 147.1, 'F': 165.2,
            'G': 75.1, 'H': 155.2, 'I': 131.2, 'K': 146.2, 'L': 131.2,
            'M': 149.2, 'N': 132.1, 'P': 115.1, 'Q': 146.2, 'R': 174.2,
            'S': 105.1, 'T': 119.1, 'V': 117.1, 'W': 204.2, 'Y': 181.2
        }
        return sum(amino_acid_weights[aa] for aa in self.sequence)

    def aminoacid_frequency(self):
        """
        Calculates frequency of aminoacids in sequence alike to nucleotide_frequency function
        """
        frequency = {}
        for amin in self.sequence:
            frequency[amin] = frequency.get(amin, 0) + 1
        return frequency

    def _validate_alphabet(self):
        valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")
        if not set(self.sequence).issubset(valid_amino_acids):
            raise ValueError("Invalid amino acid sequence")


logging.basicConfig()
logger = logging.getLogger("FastqFilter")
logger.setLevel(logging.INFO)

file_handler = logging.FileHandler("fastq_filter.log")
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)


def filter_fastq(
        input_fastq: str,
        output_fastq: str,
        gc_bounds: Union[Tuple[float, float], float] = (0, 100),
        length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
        quality_threshold: int = 0
                    ) -> Dict[str, Tuple[str, str]]:
    """
    Filters sequences in a FASTQ file based on GC content, sequence length,
    and average quality score.

    Args:
        input_fastq (str): Path to the input FASTQ file.
        output_fastq (str): Path to the output FASTQ file.
        gc_bounds (Union[Tuple[float, float], float], optional):
            GC content range for filtering. Defaults to (0, 100).
        length_bounds (Union[Tuple[int, int], int], optional):
            Length range for filtering. Defaults to (0, 2**32).
        quality_threshold (int, optional):
            Minimum average quality score. Defaults to 0.

    Returns:
        Dict[str, Tuple[str, str]]: Dictionary of filtered sequences.
    """
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)

    filtered_seqs = {}
    output_dir = 'filtered'
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, output_fastq)
    base, ext = os.path.splitext(output_path)
    counter = 1
    while os.path.exists(output_path):
        output_path = f"{base}_{counter}{ext}"
        counter += 1

    logger.info(f"File filtering started {input_fastq} → {output_path}")

    try:
        with open(output_path, 'w') as output_handle:
            for record in SeqIO.parse(input_fastq, 'fastq'):
                if not gc_filter(record.seq, gc_bounds):
                    continue
                if not length_filter(record.seq, length_bounds):
                    continue
                if not quality_filter(record.letter_annotations['phred_quality'], quality_threshold):
                    continue

                SeqIO.write(record, output_handle, 'fastq')
                filtered_seqs[record.id] = (str(record.seq), ''.join(map(chr, [q + 33 for q in record.letter_annotations['phred_quality']])))
    except Exception as e:
        logger.error(f"Error during processing {input_fastq}: {e}")
        raise

    logger.info(f"Filtering completed: {len(filtered_seqs)} selected")
    return filtered_seqs


def parse_args():
    parser = argparse.ArgumentParser(description="FASTQ filtering by GC content, length and quality")

    parser.add_argument("input_fastq", help="The path to the input FASTQ file")
    parser.add_argument("output_fastq", help="Name of the output FASTQ file")
    parser.add_argument("--gc_bounds", nargs='+', type=float, default=[0, 100],
                        help="GC content limits: one number (upper bound) or two (lower and upper bound)")
    parser.add_argument("--length_bounds", nargs='+', type=int, default=[0, 2**32],
                        help="Length limits: one number (upper bound) or two (lower and upper bound)")
    parser.add_argument("--quality_threshold", type=int, default=0,
                        help="Minimum average quality value")

    return parser.parse_args()


def quality_filter(quality_scores: list, quality_threshold: int) -> bool:
    """
    Filters based on average quality score.
    """
    if not quality_scores:
        return False
    return sum(quality_scores) / len(quality_scores) >= quality_threshold


def gc_filter(seq, gc_bounds: Tuple[float, float]) -> bool:
    """
    Filters based on GC content percentage.
    """
    gc_percent = gc_fraction(seq) * 100
    return gc_bounds[0] <= gc_percent <= gc_bounds[1]


def length_filter(seq, length_bounds: Tuple[int, int]) -> bool:
    """
    Filters based on sequence length.
    """
    return length_bounds[0] <= len(seq) <= length_bounds[1]
