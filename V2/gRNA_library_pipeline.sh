#!/bin/bash

export PATH=$PATH:/home/juan/paquetes/bowtie-1.3.1-linux-x86_64
#ejemplo, el número es la longuitud del spacer ./gRNA_library.sh glaA_ORF.fasta wrapped_glaA_CDS.fasta 20
#ejemplo ./gRNA_library_pipeline.sh Aspni_NRRL3_1_AssemblyScaffolds.fasta Aspni_NRRL3_1_GeneCatalog_CDS_20140311.fasta 18 0.7

mkdir -p gRNA_lib
# Ensure the necessary arguments are passed
if [ $# -lt 2 ]; then
    echo "Usage: $0 <ORF_FASTA> <CDS_FASTA> <SPACER_LENGTH>"
    exit 1
fi

ORF_FASTA=$1
CDS_FASTA=$2
SPACER_LENGTH=$3
CUT_OFF=$4
#GENOME_FASTA=$4


# Strip the extension from the genome file to use for directory creation
ORF_BASE=$(basename "$ORF_FASTA" .fasta)  # You can also adjust the extension (e.g., .fa)

# Create the directory for Bowtie index
mkdir -p "$ORF_BASE"_index

# Build the Bowtie index
echo "Building Bowtie index for the genome..."
bowtie-build "$ORF_FASTA" "$ORF_BASE"_index/ORF_index


# Step 2: Run Python script with input_fasta and pipe output to Bowtie
echo "Running Python script on $INPUT_FASTA and piping into Bowtie..."
python3 re_new_gRNAs_CDS.py "$CDS_FASTA" -l $SPACER_LENGTH
#bowtie alignment, -v 0 is for perfect matches, -k 100 is report multiple alignments, -f for FASTA input
echo "Searching for CDS/Genome spacer intersection through perfect alignment..."
bowtie -v 0 -k 100 -f -x "$ORF_BASE"_index/ORF_index gRNA_lib/2_all_spacers_in_CDS.fasta --al gRNA_lib/3_spacers_in_CDS_and_ORF.fasta --un gRNA_lib/4_inbetween_exons_spacers.fasta
echo "Removed CDS spacers located in exon-exon boundaries."
# Step 3: look for off-targets
echo "Allowing alignment with mismatches at the first 9 nt of the spacer..."
bowtie -n 3 --seedlen 9 -k 100 -f -x "$ORF_BASE"_index/ORF_index gRNA_lib/3_spacers_in_CDS_and_ORF.fasta -S gRNA_lib/bowtie_output.sam
awk '$0 ~ /^@/ || /XA:i:[1-9]/' gRNA_lib/bowtie_output.sam > gRNA_lib/multi_mapped.sam  #alignemnts with XA:i:1 or greater than 1 means alignments in other locations. XA:i:0 is unique alignment.

#Step 4: reject spacers with mismatched alignments (off-targets)
python3 offtargets.py gRNA_lib/multi_mapped.sam gRNA_lib/3_spacers_in_CDS_and_ORF.fasta gRNA_lib/5_no_potential_offtargets.fasta

#Step 5:
echo "Selecting a top scoring spacer for each gene Id..."
python3 filtergRNAs.py gRNA_lib/5_no_potential_offtargets.fasta gRNA_lib/6_top_scoring_spacers.fasta

#Step 6:
echo "Applying minimal score cut-off..."
python3 cutoff_filter_gRNAs.py gRNA_lib/6_top_scoring_spacers.fasta gRNA_lib/7_score_cutoff_spacers.fasta $CUT_OFF
