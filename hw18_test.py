import os
import pytest
from Bio import SeqIO
from tempfile import TemporaryDirectory
from BioinformaticToolsInator import filter_fastq

TEST_FASTQ = "data/fastq_data.fastq"

# temporary path for the FASTQ output file
@pytest.fixture
def temp_output_file():
    with TemporaryDirectory() as tmpdir:
        yield os.path.join(tmpdir, "output.fastq")


# Check: by default, filtering does not discard all reads.
def test_default_filter_does_not_discard_all(temp_output_file):
    result = filter_fastq(TEST_FASTQ, os.path.basename(temp_output_file))
    assert isinstance(result, dict)
    assert len(result) > 0, "None of the entries passed the default filtering"


# Check: too high GC threshold discards all reads
def test_gc_bounds_exclude_all(temp_output_file):
    result = filter_fastq(TEST_FASTQ, os.path.basename(temp_output_file), gc_bounds=(90, 100))
    assert len(result) == 0, "All entries would be filtered by high GC content."


# Check: only reads up to 45 nucleotides long are filtered
def test_length_bounds_include_only_short(temp_output_file):
    result = filter_fastq(TEST_FASTQ, os.path.basename(temp_output_file), length_bounds=(0, 45))
    assert all(len(seq[0]) <= 45 for seq in result.values()), "Reads longer than 45 nucleotides were found"


# Check: if the length range is too large, all reads are excluded.
def test_length_bounds_exclude_all(temp_output_file):
    result = filter_fastq(TEST_FASTQ, os.path.basename(temp_output_file), length_bounds=(200, 300))
    assert len(result) == 0, "All reads should be excluded due to the length."


# Check: filtering by a high quality threshold discards more reads
def test_quality_threshold_excludes_some(temp_output_file):
    result_low = filter_fastq(TEST_FASTQ, os.path.basename(temp_output_file), quality_threshold=0)
    result_high = filter_fastq(TEST_FASTQ, os.path.basename(temp_output_file), quality_threshold=40)
    assert len(result_high) < len(result_low), "Filtering by quality did not reduce the number of reads"


# Check: the output file is indeed created in the "filtered" folder
def test_output_file_created(temp_output_file):
    output_path = os.path.join('filtered', os.path.basename(temp_output_file))
    filter_fastq(TEST_FASTQ, os.path.basename(temp_output_file))
    assert os.path.exists(output_path), f"File {output_path} was nort created"


# Check: when specifying a single number for gc_bounds, it is perceived as the maximum
def test_single_number_gc_bound(temp_output_file):
    result = filter_fastq(TEST_FASTQ, os.path.basename(temp_output_file), gc_bounds=60)
    for seq, _ in result.values():
        gc_content = 100 * (seq.count('G') + seq.count('C')) / len(seq)
        assert gc_content <= 60, "GC content exceeds the set threshold"
