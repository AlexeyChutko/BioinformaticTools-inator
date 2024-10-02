from typing import Dict, List, Union


def transcribe(*args: str) -> Union[List[str], str]:
    results = []
    for i in args:
        transcribed_seq = i.replace('T', 'U').replace('t', 'u')
        results.append(transcribed_seq)
    if len(args) == 1:
        return "".join(results)
    return results


def reverse(*args: str) -> Union[List[str], str]:
    results = [i[::-1] for i in args]
    if len(args) == 1:
        return "".join(results)
    return results


def complement(*args: str) -> Union[List[str], str]:
    dna_dict_up = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
    dna_dict_low = {'a': 't', 'g': 'c', 't': 'a', 'c': 'g'}
    rna_dict_up = {'A': 'U', 'G': 'C', 'U': 'A', 'C': 'G'}
    rna_dict_low = {'a': 'u', 'g': 'c', 'u': 'a', 'c': 'g'}
    results = []
    for j in args:
        complement_seq = []
        if 'u' in j or 'U' in j:
            for m in j:
                if m in rna_dict_up:
                    complement_seq.append(rna_dict_up[m])
                elif m in rna_dict_low:
                    complement_seq.append(rna_dict_low[m])
                else:
                    complement_seq.append(m)
        else:
            for m in j:
                if m in dna_dict_up:
                    complement_seq.append(dna_dict_up[m])
                elif m in dna_dict_low:
                    complement_seq.append(dna_dict_low[m])
                else:
                    complement_seq.append(m)
        results.append("".join(complement_seq))
    if len(args) == 1:
        return "".join(results)
    return results


def reverse_complement(*args: str) -> Union[List[str], str]:
    normal_complement = complement(*args)
    results = normal_complement[::-1]
    if len(args) == 1:
        return "".join(results)
    return results


# Дополнительные функции
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
            results[arg] = "отсутствует"
    return results
