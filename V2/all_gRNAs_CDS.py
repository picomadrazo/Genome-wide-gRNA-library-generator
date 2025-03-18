import argparse
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq


#This step starts with a CDS catalogue and outpus all gRNA spacers that can introduce stop codons by
#cytidine deaminase base editors in the first 2-10 nt of the spacer. The first number is a score (based on
#the presence of homopolymer stretches, GC content and distance from start codon or stop codon (Sd)). The second number is Sd.

#example of usage: 
#python3 gRNAs_CDS.py wrapped_glaA_CDS.fasta 

#---------------------------------------------------------------Preparing CDS fasta file. Unwraps and adds reverse complement in next gene entry

def is_wrapped(fasta_file):
    """Check if the FASTA file contains wrapped sequences by looking for line breaks in sequence data."""
    with open(fasta_file, "r") as file:
        is_sequence_line = False
        for line in file:
            line = line.strip()
            if line.startswith(">"):  # If line starts with '>', it's a header
                is_sequence_line = False
            elif line:  # If it's a sequence line
                if is_sequence_line:  # If this is a sequence line after another sequence line
                    return True
                is_sequence_line = True
    return False

def unwrap_fasta(input_fasta, output_fasta):
    """Unwraps a FASTA file by converting multiline sequences into a single line format."""
    with open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            output_handle.write(f">{record.id}_+_\n{record.seq}\n")  # Writes header + full sequence on one line
 
def extract_reverse_complement(input_fasta, output_fasta):
    """Extract reverse complement for each sequence and add it as a new line with 'Antisense'."""
    with open(output_fasta, "w") as output_handle:  # Open in write mode to overwrite the file
        for record in SeqIO.parse(input_fasta, "fasta"):
            # Write original (sense) sequence
            output_handle.write(f">{record.id}_+_\n{record.seq}\n")
            
            # Write reverse complement (antisense) sequence
            antisense_seq = record.seq.reverse_complement()
            output_handle.write(f">{record.id}_-_\n{antisense_seq}\n")

# Argument parser
parser = argparse.ArgumentParser(description="Check and unwrap a FASTA file if wrapped, then extract reverse complement.")
parser.add_argument("input_fasta", help="Path to the input FASTA file")  # Required
parser.add_argument("-l", "--spacer_length", help="Number of nt for spacer", type=int, required=True)  # Optional flag


# Generate default output filename based on input name
args, _ = parser.parse_known_args()  # Parse only input_fasta first to get its name
input_basename = os.path.splitext(os.path.basename(args.input_fasta))[0]
default_output = f"gRNA_lib/1_{input_basename}_with_antisense.fasta"

parser.add_argument("-o", "--output", default=default_output,
                    help=f"Path to the output FASTA file (default: {default_output})")

# Parse full arguments
args = parser.parse_args()

# Check if the file is wrapped
if is_wrapped(args.input_fasta):
    print("FASTA file is wrapped. Unwrapping now...")
    unwrap_fasta(args.input_fasta, args.output)
else:
    print("FASTA file is already unwrapped. No changes made.")

# Extract reverse complement and write both sense and antisense sequences to the output file
extract_reverse_complement(args.input_fasta, args.output)
print(f"Original and reverse complement sequences added to {args.output}")



#---------------------------------------------------------------Search for spacer candidates containing codons sensitive to introduction of stop condon by C->T transitions
file = open(f"gRNA_lib/1_{input_basename}_with_antisense.fasta")
lines = []
for l in file.readlines():
    lines.append(l)

dic = {}
for i in lines:
    if re.match("^>", i):
        title = i.strip()
        dic[title] = []
    else:
        dic[title].append(i.strip("\n"))

output = open("gRNA_lib/2_all_spacers_in_CDS.fasta", "w")

for k in dic:
    sequence = "".join(dic[k])  # Join multiple lines into a single sequence
    pattern = f"(?=(.{{{args.spacer_length}}}.GG))"
    res = re.findall(pattern, sequence)

    if res:
        for idx, i in enumerate(res):  # Enumerate to create unique identifiers
            PD1 = PA1 = PB1 = PC1 = 0  # Reset position variables
            if "_-_" in k:
                In = len(sequence) - sequence.index(i)
                POS = (In - 6) / len(sequence)
                if "CCA" in i[0:11]:
                    PD = In - i[0:11].index("CCA")
                    PD1 = PD if PD % 3 == 0 else 0
            else:
                In = sequence.index(i)
                POS = (In + 17) / len(sequence)
                if "CGA" in i[1:11]:
                    PA = In + i[1:11].index("CGA") + 4
                    PA1 = PA if PA % 3 == 0 else 0
                if "CAA" in i[1:11]:
                    PB = In + i[1:11].index("CAA") + 4
                    PB1 = PB if PB % 3 == 0 else 0
                if "CAG" in i[1:11]:
                    PC = In + i[1:11].index("CAG") + 4
                    PC1 = PC if PC % 3 == 0 else 0

            # Position formula
            Sd = 2 * POS if POS <= 0.50 else 5.1535 * (POS ** 2) - 9.6508 * POS + 4.5061

            # Homopolymer penalty
            filters = ["AAAA", "TTTT", "CCCC", "GGGG"]
            shp = 0 if any(f in i for f in filters) else 1

            # GC content formula
            Sub = i[:args.spacer_length]
            GC = (Sub.count('G') + Sub.count('C')) / args.spacer_length
            Sgc = GC * 2 - 0.20 if GC < 0.60 else 2.53 - GC * 2.55

            # Final scoring
            lj = format((Sd * 0.75 + Sgc * 0.13 + shp * 0.12), '.2f')

            # Unique identifier for each sequence
            if any([PD1, PA1, PB1, PC1]):
                new_id = f"{k}i{idx}_Sc{lj}_Pos{format(POS, '.4f')}"
                output.write(f"{new_id}\n{i}\n")

output.close()

