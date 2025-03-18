import sys

def extract_qnames(sam_file):
    """Extract QNAMEs from SAM file where XA:i:1 or higher is present."""
    qnames = set()
    with open(sam_file, 'r') as sam:
        for line in sam:
            if line.startswith('@'):  # Skip header lines
                continue
            fields = line.strip().split('\t')
            qnames.add(fields[0])  # QNAME is the first column
    return qnames

def filter_fasta(fasta_file, qnames, output_fasta):
    """Remove entries from FASTA file that have QNAMEs in the set."""
    with open(fasta_file, 'r') as fasta, open(output_fasta, 'w') as out_fasta:
        write_seq = True
        for line in fasta:
            if line.startswith(">"):
                qname = line.split()[0][1:]  # Extract QNAME from FASTA header
                write_seq = qname not in qnames  # Check if QNAME is in the exclusion list
            if write_seq:
                out_fasta.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <multi_mapped.sam> <input.fasta> <filtered_output.fasta>")
        sys.exit(1)

    sam_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_fasta = sys.argv[3]

    qnames_to_remove = extract_qnames(sam_file)
    filter_fasta(fasta_file, qnames_to_remove, output_fasta)

    print(f"Filtered FASTA saved as: {output_fasta}")
    print(f"Removed {len(qnames_to_remove)} gene entries with potential off-targets.")
    
#usage python3 offtargets.py gRNA_lib/multi_mapped.sam gRNA_lib/3_spacers_in_CDS_and_ORF.fasta gRNA_lib/5_no_potential_offtargets_.fasta

