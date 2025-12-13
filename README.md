# Mutation Rate and Recombination Analysis Pipeline
A reproducible bioinformatics scripts for estimating mutation rates and analyzing recombination and GC content correlations in *Bos taurus* and *Bos indicus hybrid*, using Ensembl whole-genome multiple alignments (EPO).
This pipeline retrieves multiple whole-genome alignments (MAF), extracts mutations, filters for putatively neutral sites, calculates relative and absolute mutation rates, and examines recombination and GC content correlations. 
All analyses use open-source tools such as Biopython, BEDTools, and R.

##  Project Structure

The repository is organized as follows:

genome-mutation-pipeline/
│
├── README.md                # ← Main project description and instructions
│
├── data/                    # Raw and processed data 
│   ├── alignment.maf
│   └── gene_annotations/
│
├── scripts/                 # All Python and R scripts 
│   ├── 01_filter_alignments.py
│   ├── 02_extract_mutations.py
│   ├── 03_filter_neutral_sites.sh
│   ├── 04_count_mutations.py
│   ├── 05_count_substitution_mutations.py
│   ├── 06_extract_trinucleotide_contexts
│   ├── 07_recombination_rate_estimation
│   └── 08_GC_content.sh
│
├── results/                 
│   ├── mutation_rates.csv
│   ├── plots/


## Installation
Clone this repository and install dependencies:
```bash
# Download your project
git clone https://github.com/Redempta-Kajungiro/MSc-Bioinformatics-2025_cattle-mutation-rate-pipeline.git
cd genome-mutation-pipeline

# Create environment
conda env create -f environment.yml

# Activate environment
conda activate genome-mutation
## Usage
1. Preprocess alignments  
   ```bash
   python scripts/01_filter_alignments.py

2. Extract mutations
   scripts/02_extract_mutations.py
   
3. Filter neutral sites
   scripts/03_filter_neutral_sites.sh

4. Calculate mutation rates

    scripts/04_count_mutations.py
            05_count_substitution_mutations.py
6.Tri-nucleotide context matrix
      scripts/06_extract_trinucleotide_contexts.py
7. Recombination and GC content
      scripts/07_recombination_rate_estimation.py
      scripts/08_GC_content.sh

## Outputs
- `results/tables/` → summary tables (mutation and recombination)
- `results/plots/` → heatmaps, scatterplots
- `results/logs/` → run logs
## Reproducibility and Adaptability
The pipeline is modular and can be adapted to any species with available Ensembl EPO multiple alignments.
It relies entirely on open-source tools (Biopython, BEDTools, R) and standard file formats (MAF, BED, GTF), ensuring reproducibility across environments.
Scripts can be customized for new genome assemblies by updating the input MAF and GTF files.
