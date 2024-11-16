from typing import Dict, List, Union

dna_dict = {
    'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G',
    'a': 't', 'g': 'c', 't': 'a', 'c': 'g'
    }
rna_dict = {
    'A': 'U', 'G': 'C', 'U': 'A', 'C': 'G',
    'a': 'u', 'g': 'c', 'u': 'a', 'c': 'g'
    }


def is_valid_na(seq: str) -> None:
    if ('u' in seq or 'U' in seq) and ('t' in seq or 'T' in seq):
        raise ValueError("Error: simultaneous presence of 'u' and 't' "
                         " in any register is unacceptable.")
    if not set([symb.upper() for symb in seq]).issubset(
            {'A', 'T', 'G', 'C', 'U'}
            ):
        raise ValueError(f"Error: the sequence contains "
                         f"invalid characters: {seq}")


def transcribe(*args: str) -> Union[List[str], str]:
    results = []
    for seq in args:
        transcribed_seq = seq.replace('T', 'U').replace('t', 'u')
        results.append(transcribed_seq)
    if len(args) == 1:
        return results[0]
    return results


def reverse(*args: str) -> Union[List[str], str]:
    results = [seq[::-1] for seq in args]
    if len(args) == 1:
        return results[0]
    return results


def complement(*args: str) -> Union[List[str], str]:
    results = []
    for seq in args:
        complement_seq = []
        if 'u' in seq or 'U' in seq:
            for symb in seq:
                if symb in rna_dict:
                    complement_seq.append(rna_dict[symb])
                elif symb in rna_dict:
                    complement_seq.append(rna_dict[symb])
                else:
                    complement_seq.append(symb)
        else:
            for symb in seq:
                if symb in dna_dict:
                    complement_seq.append(dna_dict[symb])
                elif symb in dna_dict:
                    complement_seq.append(dna_dict[symb])
                else:
                    complement_seq.append(symb)
        results.append("".join(complement_seq))
    if len(args) == 1:
        return results[0]
    return results


def reverse_complement(*args: str) -> Union[List[str], str]:
    return reverse(complement(*args))


def nucleotide_frequency(*args: str) -> Dict[str, Dict[str, int]]:
    frequency_dict = {}
    for sequence in args:
        seq = sequence.upper()
        frequency = {}
        for nucl in seq:
            if nucl in frequency:
                frequency[nucl] += 1
            else:
                frequency[nucl] = 1
        frequency_dict[sequence] = frequency
    return frequency_dict


def gc_content(*args: str) -> Dict[str, float]:
    gc_dict = {}
    for seq in args:
        seq = seq.upper()
        g_count = seq.count('G')
        c_count = seq.count('C')
        gc_count = g_count + c_count
        percent = (gc_count / len(seq)) * 100 if seq else 0
        gc_dict[seq] = int(round(percent, 0))
    return gc_dict


def find_start_codons(*args: str) -> Dict[str, Union[List, str]]:
    results = {}
    for arg in args:
        seq = arg.upper()
        start_codon = "AUG" if 'U' in seq else "ATG"
        positions = []
        for i in range(len(seq) - 2):
            if seq[i:i+3] == start_codon:
                positions.append(i)
        if positions:
            results[arg] = positions
        else:
            results[arg] = []
    return results
