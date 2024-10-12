import os


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None) -> None:
    """
    Converts a multi-line FASTA file to a one-line format for each sequence
     Args:
        input_fasta (str): The path to the input FASTA file with sequences
        output_fasta (str, optional): The path to the output FASTA file where each sequence will be written 
                                      on one line. If there is no output_fasta, the output file will be named as 
                                      the input file with '_oneline' appended before the extension.

    Returns:
        None: The converted FASTA file is saved to disk
    """
    if output_fasta is None:
        base, ext = os.path.splitext(input_fasta)
        output_fasta = f"{base}_oneline{ext}"
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        sequence = []
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    outfile.write(''.join(sequence) + '\n')
                    sequence = []
                outfile.write(line + '\n')
            else:
                sequence.append(line)
        if sequence:
            outfile.write(''.join(sequence) + '\n')
    print(f"Converted FASTA file saved to: {output_fasta}")


def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Parses BLAST output to extract the best match protein names and saves them to a file.

    Args:
        input_file (str): The path to the BLAST output file
        output_file (str): The path to the output file where the sorted protein names will be saved

    Returns:
        None: The parsed results are saved to the specified output file
    """
    best_matches = set()
    with open(input_file, 'r') as infile:
        in_significant_alignments = False
        for line in infile:
            if line.startswith("Sequences producing significant alignments:"):
                in_significant_alignments = True
                continue
            if in_significant_alignments and line.strip() == "":
                in_significant_alignments = False
                continue
            if in_significant_alignments:
                parts = line.split()
                if len(parts) > 0:
                    description = " ".join(parts[:-5])
                    if "[" in description:
                        protein_name = description.split(']')[0] + "]"
                        best_matches.add(protein_name)
    sorted_matches = sorted(best_matches)
    with open(output_file, 'w') as outfile:
        for match in sorted_matches:
            outfile.write(match + '\n')

    print(f"Parsing results saved in {output_file}")
