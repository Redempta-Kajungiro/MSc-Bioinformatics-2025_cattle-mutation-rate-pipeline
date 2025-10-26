# Mutation Rate and Recombination Analysis Pipeline
A reproducible bioinformatics workflow for estimating mutation rates and analyzing recombination and GC content correlations in *Bos taurus* and *Bos indicus hybrid*, using Ensembl whole-genome multiple alignments (EPO).
This pipeline retrieves multiple whole-genome alignments (MAF), extracts lineage-specific mutations, filters for putatively neutral sites, calculates relative and absolute mutation rates, and examines recombination and GC content correlations. All analyses use open-source tools such as Biopython, BEDTools, and R.
## Repository Structure
genome-mutation-pipeline/
├── data/
├── scripts/
├── notebooks/
├── results/
└── utils/
