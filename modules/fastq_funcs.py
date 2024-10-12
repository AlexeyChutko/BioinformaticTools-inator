from typing import Tuple, Union
import os


def read_fastq(input_fastq: str) -> dict:
    sequences = {}
    with open(input_fastq, 'r') as file:
        while True:
            name = file.readline().strip()
            if not name:
                break
            sequence = file.readline().strip()
            file.readline()
            quality = file.readline().strip()
            sequences[name] = (sequence, quality)
    return sequences


def write_fastq(sequences: dict, output_fastq: str):
    base, ext = os.path.splitext(output_fastq)
    counter = 1
    while os.path.exists(output_fastq):
        output_fastq = f"{base}_{counter}{ext}"
        counter += 1
    output_dir = 'filtered'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    with open(os.path.join(output_dir, output_fastq), mode='w') as file:
        for name, (sequence, quality) in sequences.items():
            file.write(f"{name}\n{sequence}\n+\n{quality}\n")


def quality_filter(
        quality_string: str,
        quality_threshold: int
        ) -> float:
    qual_scores = [ord(i) - 33 for i in quality_string]
    if len(qual_scores) == 0:
        mean_qual = 0
    else:
        mean_qual = sum(qual_scores) / len(qual_scores)
    return mean_qual >= quality_threshold


def gc_filter(
        seq: str,
        gc_bounds: Union[Tuple[float, float], float]
        ) -> float:
    g_count = seq.count('G') + seq.count('g')
    c_count = seq.count('C') + seq.count('c')
    gc_count = g_count + c_count
    if len(seq) == 0:
        gc_percent = 0
    else:
        gc_percent = (gc_count / len(seq)) * 100
    return gc_bounds[0] <= gc_percent <= gc_bounds[1]


def length_filter(
        seq: str,
        length_bounds: Union[Tuple[int, int], int]
        ) -> int:
    seq_length = len(seq)
    return length_bounds[0] <= seq_length <= length_bounds[1]
