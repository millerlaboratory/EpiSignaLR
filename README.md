# EpiSignaLR

## v1.0 


A workflow for methylation classification of mendelian conditions using methylation aware sequencing.

EpiSignaLR supports Oxford Nanopore (ONT), PacBio, and Illumina 5-base data as input. When running `EpiSignaLR.sh`, select the platform with `--platform`. Coordinates are for the hg38 reference genome.

This will generate a pileup, filter the pileup to the target CpGs, then run:
- Array-derived SVM classifiers for each of 34 conditions
    (Adapted from: NSBEpi, https://github.com/JorisVermeeschLab/NSBEpi/tree/main)
- LRS-derived ElasticNet+SVM classifiers for **SETD1B-related disorder**, **Kleefstra syndrome**, and **Kabuki syndrome**


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
conda activate episignalr_env

chmod +x EpiSignaLR.sh

```

## Basic Usage

```
./EpiSignaLR.sh --platform <ont|pacbio|illumina5base> --sample <SAMPLE> --outdir <OUTPUT_DIR> \
    --input <BAM|gVCF> [--ref <hg38.fa>] [--threads N] [--models SETD1B,Kleefstra,Kabuki]
```

| Option | Description |
|---|---|
| `-p, --platform` | `ont`, `pacbio`, or `illumina5base` |
| `-s, --sample` | Sample ID; used as the output prefix and subdirectory name |
| `-o, --outdir` | Output directory; results go to `<outdir>/<sample>/` |
| `-i, --input` | ONT modBAM, PacBio aligned BAM, or Illumina 5-base methylation gVCF (all hg38) |
| `--pileup` | Use an existing pileup instead of `--input` and skip the pileup step |
| `-r, --ref` | hg38 reference FASTA (required for `ont` and `pacbio` with `--input`) |
| `-t, --threads` | Threads for the pileup step (`modkit` / `aligned_bam_to_cpg_scores`); default 4 |
| `-m, --models` | Comma-separated LRS-derived models to run; default `SETD1B,Kleefstra,Kabuki` |

Examples:
```
# ONT
./EpiSignaLR.sh -p ont -s SAMPLE1 -o results/ -i SAMPLE1.bam -r hg38.fa -t 16

# PacBio (requires pb-CpG-tools' aligned_bam_to_cpg_scores on PATH)
./EpiSignaLR.sh -p pacbio -s SAMPLE2 -o results/ -i SAMPLE2.bam -r hg38.fa -t 16

# Illumina 5-base
./EpiSignaLR.sh -p illumina5base -s SAMPLE3 -o results/ -i SAMPLE3.methyl.gvcf.gz

# Existing modkit pileup
./EpiSignaLR.sh -p ont -s SAMPLE1 -o results/ --pileup SAMPLE1.modkit.pileup_combined.bed.gz
```
A directory with the sample ID is made in the given output diretory.

Steps whose outputs already exist are skipped, so an interrupted run can be resumed by rerunning the same command.
*If a run is interrupted and output files are incomplete, the interrupted step may not re-run without deleting certain output files*

### Code layout
```
EpiSignaLR.sh              # Argument parsing and orchestration
lib/common.sh              # Logging and helper functions
lib/steps.sh               # Platform-independent classification steps
lib/platforms/<name>.sh    # Per-platform pileup creation and format conversion
```

### Data
The data directory includes the methylation fractions calculated per training sample for each nanopore-sequencing derived condition classifier as well as the classifer model and metadata. 


## Steps

### 1. Pileup

**ONT**: modkit v0.6.1 is used to generate a bedMethyl pileup using the `--combine-mods` flag to handle 5hmC probaililites

**PacBio**: pb-cpg-tools v3.0.0 `aligned_bam_to_cpg_scores` with the flags `--pileup-mode count` and `--modsites-mode reference`

**Illumina5base**: CpG fraction values from the the hard-filtered.gvcf output from the DRAGEN pipeline are extracted and formatted into a bed file for efficient use in down stream steps


### 2. Array-dervied classifier scores

The pileup is intersected with CpG positions for running the array-derived classifers and formatted to match the expected output for `NSBEpi`. In addition to an output with the scores and final conclusion, a barplot of the top 5 scores is generated.

### 3. Nanopore sequencing-derived classifier scores

A weighted average methylation fraction per DMR in the condition is calculated from the pileup to be used for each condition classifer. A PCA plot of the sample projected onto the training data and a heatmap of the methylation values for the sample and training cohort are included as outlputs for visual analysis and interpretation.


## Output

The pipeline generates a sample-specific directory containing:
```
<sample_name>/
├── Pileups/
│   ├── <sample>_modkit.pileup_combined.bed.gz      # Un-phased modkit pileup, used for filtering for classification
│   ├── <sample>_modkit.pileup_hp1.bed.gz           # Phased modkit pileup for hp1, not used
│   ├── <sample>_modkit.pileup_hp1.bed.gz           # Phased modkit pileup for hp2, not used
│   └── <sample>.log                                # Modkit pileup log
├── Array_Derived_Episignature_Results/
│   ├──<sample>_modkit.episignature.bed            # Subsetted pileup to episignature positions
│   ├──<sample>_results_score.csv                  # The scores for each condition
│   ├──<sample>_results_top_result.txt             # The top score and associated condition (or control)
│   ├──<sample>_results_top5_barplot.png           # A barplot of the top 5 scores for the sample
├── LRS_Derived_Episignature_Results/
│   ├── <sample>.SETD1B_enet_svm_full_panel_DMPs.csv          # Region-averaged methylation at model DMRs
│   ├── <sample>.SETD1B_enet_svm_full_panel_DMPs_summary.txt  # Summary of regeion-average calculation
│   ├── <sample>.SETD1B.svm_results_score.csv                 # SVM score, predicted class, indeterminate flag
│   ├── <sample>.SETD1B.svm_results_top_result.txt            # Human-readable top result
│   ├── <sample>.SETD1B.svm_results_pca.svg                   # PCA: training cohort + this sample
│   └── <sample>.SETD1B.svm_results_heatmap.svg               # Heatmap: training cohort + this sample at selected features
│   ├── <sample>.Kabuki_enet_svm_full_panel_DMPs.csv          # Region-averaged methylation at model DMRs
│   ├── <sample>.Kabuki_enet_svm_full_panel_DMPs_summary.txt  # Summary of regeion-average calculation
│   ├── <sample>.Kabuki.svm_results_score.csv                 # SVM score, predicted class, indeterminate flag
│   ├── <sample>.Kabuki.svm_results_top_result.txt            # Human-readable top result
│   ├── <sample>.Kabuki.svm_results_pca.svg                   # PCA: training cohort + this sample
│   └── <sample>.Kabuki.svm_results_heatmap.svg               # Heatmap: training cohort + this sample at selected features
│   ├── <sample>.Kleefstra_enet_svm_full_panel_DMPs.csv          # Region-averaged methylation at model DMRs
│   ├── <sample>.Kleefstra_enet_svm_full_panel_DMPs_summary.txt  # Summary of regeion-average calculation
│   ├── <sample>.Kleefstra.svm_results_score.csv                 # SVM score, predicted class, indeterminate flag
│   ├── <sample>.Kleefstra.svm_results_top_result.txt            # Human-readable top result
│   ├── <sample>.Kleefstra.svm_results_pca.svg                   # PCA: training cohort + this sample
│   └── <sample>.Kleefstra.svm_results_heatmap.svg               # Heatmap: training cohort + this sample at selected features
```

*EpiSignaLR is intended for research use ONLY, not for diagnostic procedures*