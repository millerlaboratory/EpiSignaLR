# EpiSignaLR

A workflow for methylation classification on 34 mendelian conditions. Adapted from: https://github.com/JorisVermeeschLab/NSBEpi/tree/main

To run, use the EpiSignaLR.sh script which requires a BAM, a sample ID ($PREFIX) and an output directory.

This will generate a pileup, filter the pileup to the target CpGs, and then run SVM classifiers for each condition

The outputs will be the full pileup, the filtered pileup, all 35 scores, the top score in txt file, and a png barplot of the top 5 scores