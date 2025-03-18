# Genome-wide-gRNA-library generator
Tool to generate genome-wide gRNA libraries for fungal genomes.

The starting material are genomic sequences (genome assembly & CDS catalogues) and the output is a .fasta file with a single spacer for each predicted gene. 

Usage is
./gRNA_library_pipeline.sh <GENOME_FASTA> <CDS_FASTA> <SPACER_LENGTH> <CUT_OFF>
  <GENOME_FASTA>  Genome assmebly file.
  <CDS_FASTA>  Coding DNA sequence file, i.e. no introns.
  <SPACER_LENGTH>  Desired spacer length. Recommended: 18 nt.
  <CUT_OFF> From 0.0 to 1.0 Based on internal formula considering distance from the start codon, GC content, and presence of homopolymer tracts. Recommended: 0.7

Gather the. fasta files and bash script in a dedicated directory and run. Requires Bowtie 1 available in PATH. 
Example: ./gRNA_library_pipeline.sh Aspni_NRRL3_1_AssemblyScaffolds.fasta Aspni_NRRL3_1_GeneCatalog_CDS_20140311.fasta 18 0.7

Working example is _Aspergillus niger_ NRRL3.


Created by Jun Lyu and Juan Pablo Morán Torres for Utrecht University Microbiology lab under the supervision of H. A. B. Wösten.
