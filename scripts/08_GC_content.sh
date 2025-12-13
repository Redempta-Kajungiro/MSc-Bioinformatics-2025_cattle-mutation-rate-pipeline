#!/bin/bash
# Required software
#   - samtools
#   - bedtools

# conda install -c bioconda samtools bedtools

# Step 1: Index the genome
#========================
GENOME_FASTA="genome.fasta"   # <- Replace with your genome file
samtools faidx $GENOME_FASTA
# This creates genome.fasta.fai

# Step 2: Create genome size file

GENOME_SIZE="genome.genome"
cut -f1,2 ${GENOME_FASTA}.fai > $GENOME_SIZE
# This file contains chromosome names and lengths

# Step 3: Create fixed-size bins
BIN_SIZE=100000  # 100 kb bins (adjust as needed)
BINS_FILE="genome_bins_${BIN_SIZE}bp.bed"
bedtools makewindows -g $GENOME_SIZE -w $BIN_SIZE > $BINS_FILE

# Step 4: Remove bins overlapping genes

GENES_BED="genes.bed"  # Replace with your BED file of genes
NEUTRAL_BINS="neutral_bins_${BIN_SIZE}bp.bed"
bedtools intersect -v -a $BINS_FILE -b $GENES_BED > $NEUTRAL_BINS

# Step 5: Calculate GC content per bin
GC_OUTPUT="gc_content_${BIN_SIZE}bp.txt"
bedtools nuc -fi $GENOME_FASTA -bed $NEUTRAL_BINS > $GC_OUTPUT

# Step 6: Count mutations per bin
MUTATIONS_BED="filtered_mutations.bed"  # Replace with your mutation BED file
MUTATION_COUNTS="mutation_counts_${BIN_SIZE}bp.tsv"
bedtools intersect -a $NEUTRAL_BINS -b $MUTATIONS_BED -c > $MUTATION_COUNTS

echo "Workflow complete!"
echo "Bins: $NEUTRAL_BINS"
echo "GC content: $GC_OUTPUT"
echo "Mutation counts: $MUTATION_COUNTS"

