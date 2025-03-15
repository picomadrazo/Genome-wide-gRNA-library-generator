import argparse
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq

#---------------------------------------------------------------Preparing CDS fasta file. Unwraps

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
            output_handle.write(f">{record.id}\n{record.seq}\n")  # Writes header + full sequence on one line
 
def extract_reverse_complement(input_fasta, output_fasta):
    """Extract reverse complement for each sequence and add it as a new line with 'Antisense'."""
    with open(output_fasta, "w") as output_handle:  # Open in write mode to overwrite the file
        for record in SeqIO.parse(input_fasta, "fasta"):
            # Write original (sense) sequence
            output_handle.write(f">{record.id}\n{record.seq}\n")
            
            # Write reverse complement (antisense) sequence
            antisense_seq = record.seq.reverse_complement()
            output_handle.write(f">{record.id}_Antisense\n{antisense_seq}\n")

# Argument parser
parser = argparse.ArgumentParser(description="Check and unwrap a FASTA file if wrapped, then extract reverse complement.")
parser.add_argument("input_fasta", help="Path to the input FASTA file")  # Required

# Generate default output filename based on input name
args, _ = parser.parse_known_args()  # Parse only input_fasta first to get its name
input_basename = os.path.splitext(os.path.basename(args.input_fasta))[0]
default_output = f"{input_basename}_with_antisense.fasta"

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
file = open(f"{input_basename}_with_antisense.fasta")
lines =[]
for l in file.readlines():
    lines.append(l)

dic={}
for i in lines:
    if re.match("^>",i):
        title = i.strip()
        dic[title]=""
    else:
        dic[title] = dic[title] + i.strip("\n")

fif = []
output=open("Result heading scoring with distance.txt","w")  #last number is the distance from the start codon (or from stop codon in antisense)
for k in dic:
    res=re.findall("(?=(.{20}.GG))",dic[k])
    if res:
        output.write("\n" + k+"\n")
        for i in res:
            PD1 = 0
            PA1 = 0
            PB1 = 0
            PC1 = 0
            if "Antisense" in k:
                In = len(dic[k]) - dic[k].index(i)
                POS = (In - 6)/len(dic[k])
                if "CCA" in i[0:11]:
                    PD = In - i[0:11].index("CCA")
                    if PD % 3 == 0:
                        PD1 = PD
                    else:
                        PD1 = 0
                else:
                    PD1 = 0
            else:
                In = dic[k].index(i)
                POS = (In + 17)/len(dic[k])
                if "CGA" in i[1:11]:
                    PA = In + i[1:11].index("CGA") + 4
                    if PA % 3 == 0:
                        PA1 = PA
                    else:
                        PA1 = 0
                else:
                    PA1 = 0
                if "CAA" in i[1:11]:
                    PB = In + i[1:11].index("CAA") + 4
                    if PB % 3 == 0:
                        PB1 = PB
                    else:
                        PB1 = 0
                else:
                    PB1 = 0
                if "CAG" in i[1:11]:
                    PC = In + i[1:11].index("CAG") + 4
                    if PC % 3 == 0:
                        PC1 = PC
                    else:
                        PC1 = 0
                else:
                    PC1 = 0
            if POS <= 0.50:  #position formula
                Sd = 2 * POS
            else:
                Sd= 5.1535*(POS*POS)-(9.6508*POS)+4.5061 #Sd = 2-2 * POS
            filters = ["AAAA","TTTT","CCCC","GGGG"] #homopolymer penalty component
            if any(name in i for name in filters) == False:
                shp = 1
            else:
                shp = 0
            Sub = i[:20]
            GC = (Sub.count('G') + Sub.count('C'))/ 20
            if GC < 0.60: #GC content formula
                Sgc = GC * 2 - 0.20
            else:
                Sgc = 2.53 - GC * 2.55
            lj = format((Sd * 0.75 + Sgc * 0.13 + shp * 0.12), '.2f')
            if PD1 !=0 or PA1 != 0 or PB1 != 0 or PC1 != 0 or PD1 != 0:
                fif = str(lj) + "," +  i + ","+ format((POS),'.2f') + "\n"
                output.write(str(fif))
output.close()
