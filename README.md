# Mutation Rate and Recombination Analysis Pipeline
A reproducible bioinformatics workflow for estimating mutation rates and analyzing recombination and GC content correlations in *Bos taurus* and *Bos indicus hybrid*, using Ensembl whole-genome multiple alignments (EPO).
This pipeline retrieves multiple whole-genome alignments (MAF), extracts lineage-specific mutations, filters for putatively neutral sites, calculates relative and absolute mutation rates, and examines recombination and GC content correlations. All analyses use open-source tools such as Biopython, BEDTools, and R.

##  Project Structure

The repository is organized as follows:

genome-mutation-pipeline/
│
├── environment.yml          # ← The Conda environment setup file 
├── README.md                # ← Main project description and instructions
├── LICENSE                  # ← (optional but recommended)
│
├── data/                    # Raw and processed data 
│   ├── example_alignment.maf
│   └── gene_annotations/
│
├── scripts/                 # All Python and R scripts for workflow
│   ├── 01_preprocess_maf.py
│   ├── 02_extract_mutations.py
│   ├── 03_filter_neutral_sites.py
│   ├── 04_calculate_mutation_rate.R
│   └── 05_gc_recombination_analysis.py
│
├── results/                 # Output files (plots, tables, etc.)
│   ├── mutation_rates.csv
│   ├── plots/
│   └── logs/
│
└── workflow/                # (optional) Workflow diagram, Snakemake, or Nextflow scripts
    └── workflow_diagram.png

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
   python scripts/02_filter_alignments.py

2. Extract mutations

    python scripts/03_extract_mutations.py
3. Filter neutral sites
  scripts/05_filter_neutral_sites.py

4. Calculate mutation rates

   python scripts/06_mutation_rate_estimation.py

5. Visualize results in R
   Open notebooks/mutation_rate_visualization.Rmd
#### 📊 6. **Outputs**
List what users will find after running the workflow.
```markdown
## Outputs
- `results/tables/` → summary tables (mutation and recombination)
- `results/plots/` → heatmaps, scatterplots
- `results/logs/` → run logs
## Reproducibility and Adaptability
The pipeline is modular and can be adapted to any species with available Ensembl EPO multiple alignments.
It relies entirely on open-source tools (Biopython, BEDTools, R) and standard file formats (MAF, BED, GTF), ensuring reproducibility across environments.
Scripts can be customized for new genome assemblies by updating the input MAF and GTF files.
