from typing import Tuple, Union


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
