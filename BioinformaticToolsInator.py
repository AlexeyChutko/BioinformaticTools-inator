from typing import Dict, List, Tuple, Union
from modules import dna_rna_funcs, fastq_funcs


def run_dna_rna_tools(
        *args: str
        ) -> Union[List[str], Dict[str, Dict[str, int]], str]:
    """
    Performs a specified action on all sequences transmitted from input

    Args:
    - DNA or RNA sequences as str and a procedure to perform
      as the last argument

    Possible procedures:
    - Main procedures: reverse, transcribe, complement, reverse_complement
    - Additional procedures:
      nucleotide_frequency, gc_content, find_start_codons

    Returns:
    - Union[List[str], Dict[str, Dict[str, int]], str]:
      sequence or a collection of sequences with paramters
    """
    for j in args[:-1]:
        if ('u' in j or 'U' in j) and ('t' in j or 'T' in j):
            raise ValueError("Error: simultaneous presence of 'u' and 't' "
                             " in any register is unacceptable.")
        if not set([i.upper() for i in j]).issubset({'A', 'T', 'G', 'C', 'U'}):
            raise ValueError(f"Error: the sequence contains "
                             f"invalid characters: {j}")
    if args[-1] == 'reverse':
        return dna_rna_funcs.reverse(*args[:-1])
    elif args[-1] == 'transcribe':
        return dna_rna_funcs.transcribe(*args[:-1])
    elif args[-1] == 'complement':
        return dna_rna_funcs.complement(*args[:-1])
    elif args[-1] == 'reverse_complement':
        return dna_rna_funcs.reverse_complement(*args[:-1])
    elif args[-1] == 'nucleotide_frequency':
        return dna_rna_funcs.nucleotide_frequency(*args[:-1])
    elif args[-1] == 'gc_content':
        return dna_rna_funcs.gc_content(*args[:-1])
    elif args[-1] == 'find_start_codons':
        return dna_rna_funcs.find_start_codons(*args[:-1])
    else:
        print('Unidentified instruction')


def filter_fastq(
        seqs: Dict[str, Tuple[str, str]],
        gc_bounds: Union[Tuple[float, float], float] = (0, 100),
        length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
        quality_threshold: int = 0
        ) -> Dict[str, Tuple[str, str]]:
    """
    Filters provided sequences based on its gc composition, length
    and average sequence quality

    Args:
    - seqs Dict[str, Tuple[str, str]]: a dictionary consisting of
      fastq sequents.
        The key is a string, the name of the sequence.
        The value is a tuple of two strings: sequence and quality.
    - gc_bounds Union[Tuple[float, float], float]:
      the GC interval of the composition (in percent) for filtration
    - length_bounds Union[Tuple[int, int], int]:
      the length interval for filtering
      within which the filtered sequences should be
    - quality_threshold int:
      the threshold value of the average read quality for filtering

    Returns:
    - Dict[str, Tuple[str, str]]: a dictionary with filtered sequences
    """

    filtered_seqs = {}
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)

    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)

    for name, (sequence, quality) in seqs.items():
        if not fastq_funcs.gc_filter(sequence, gc_bounds):
            continue
        if not fastq_funcs.length_filter(sequence, length_bounds):
            continue
        if not fastq_funcs.quality_filter(quality, quality_threshold):
            continue
        filtered_seqs[name] = (sequence, quality)

    return filtered_seqs
