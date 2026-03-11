# EpiSignaLR

A workflow for methylation classification on 35 mendelian conditions. Adapted from: https://github.com/JorisVermeeschLab/NSBEpi/tree/main

To run, use the EpiSignaLR.sh script which requires a BAM, a sample ID ($PREFIX), the output directory, and a reference genome. Coordinates are for the hg38 reference genome.

This will generate a pileup, filter the pileup to the target CpGs, and then run SVM classifiers for each condition

The outputs will be the full pileup, the filtered pileup, all 35 scores, the top score in txt file, and a png barplot of the top 5 scores

## Installation

### Prerequisites

-Linux/Unix environment

-Conda or Mamba package manager

-R (>= 4.0)

-Python (>= 3.13)

### Quick Install

```
bash

# Clone the repository

git clone git@github.com:millerlaboratory/EpiSignaLR.git
cd EpiSignaLR


# Create conda environment

conda env create -f environment.yml
conda activate episignalr-env

chmod +x EpiSignaLR.sh

```

# Basic Usage
**EpisignaLR takes a bam file aligned to hg38 and an hg38 reference genome to generate a methylation pileup and then classify based on the mathylation fractions**
```
bash ./TRoLR.sh <BAM> [SAMPLE NAME] [OUTPUT_DIR] <REF>
```

## Output

The pipeline generates a sample-specific directory containing:
```
<sample_name>/
├── <sample>_modkit.pileup.bed.gz      # Modkit pileup
├── <sample>_modkit.episignature.bed   #Subsetted pileup to episignature positions
├── <sample>_results_score.csv         # The scores for each condition
├── <sample>_results_top_result.txt    # The top score and associated condition (or control)
├── <sample>_results_top5_barplot.png   # A barplot of the top 5 scores for the sample
```