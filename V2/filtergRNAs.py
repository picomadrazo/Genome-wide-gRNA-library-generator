import re
import sys
from Bio import SeqIO

def extract_gene_id(header):
    """Extract the gene ID while ignoring strand direction and other fields."""  #>jgi|Aspni_NRRL3_1|32|NRRL3_00032_+_i36_Sc0.80_Pos0.5531
    match = re.search(r"(NRRL3_\d+_)", header)
    return match.group(1) if match else None

def extract_score(header):
    """Extract the Sc score from the FASTA header."""
    match = re.search(r"Sc([\d.]+)", header)
    return float(match.group(1)) if match else None

def filter_highest_sc(fasta_file, output_file):
    best_sequences = {}

    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            gene_id = extract_gene_id(record.id)
            score = extract_score(record.id)

            if gene_id and score is not None:
                if gene_id not in best_sequences or score > best_sequences[gene_id][0]:
                    best_sequences[gene_id] = (score, record)

    with open(output_file, "w") as output_handle:
        for score, record in best_sequences.values():
            SeqIO.write(record, output_handle, "fasta")

# Usage
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python filtergRNAs.py input.fasta output.fasta")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]

    filter_highest_sc(input_fasta, output_fasta)
    print(f"Filtered sequences saved to {output_fasta}")

