# TACIT (Target Chromatin Indexing and Tagmentation)

## TACIT analysis
tacit.sh - a script for analyzing TACIT data incuding mapping to mm10 reference genome and removing PCR duplicates.

## CoTACIT analysis
1. Split raw sequening data for each histone modification based on T5 and T7 barcodes, with scripts located in the CoTACIT folder.
2. Analyze sequencing data for each histone modification using tacit.sh.

## scRNA-seq analysis
1. Split raw sequening data for each single cell based on oligo dT barcodes, with scripts in the scRNA/00_split single cells folder.
2. Analyze sequencing data including cutting adaptor, mapping to mm10 reference genome, and calculating TPM, with scripts in the scRNA/01_scRNA-QC folder.

## ChromHMM analysis
1. ChromHMM manual: https://ernstlab.biolchem.ucla.edu/ChromHMM/
2. To reduce noise and mitochondrial interference, all reads from mitochondrial chromatin are filtered out.
3. For Fig.4, all ChromHMM model and emission probabilities txt file for synthetic cells are listed in chromhmm/toti folder. And the genome-wide chromatin states annotation at a 200-kb resolution for all cells are provided.
4. For Fig.5, the forward-backward algorithm is used to learn the posterior probability distribution for all synthetic cells, with a model txt file provided in the chromhmm/ICM_TE folder. And the genome-wide chromatin states annotation at a 2000-kb resolution for all cells are also provided.

## Integration analysis
#### 1. Multimodality integration with scRNA-seq as anchors follows the Seurat V4 anchoring workflow: https://satijalab.org/seurat/.
#### 2. WNN analysis: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.
#### 3. scDeepHIM:
(1) Read count matrices for histone modifications for each samples in training groups were generated using either the bamCoverage function or an alternative method. These matrices were partitioned into input files (H3K27ac, H3K27me3, H3K9me3) and target files (H3K4me1, H3K4me3, H3K36me3). We trained scDeepHIM models with profiles of six histone modifications using 80% of the samples and evaluated its performance on the remaining 20%. During the training phase, the scDeepHIM model optimized an objective function to predict H3K4me1, H3K4me3, and H3K36me3 signals using the provided H3K27ac, H3K27me3, and H3K9me3 signals from both bulk ChIP-seq profiles and single-cell ChIP-seq profiles.
(2) Following training, scDeepHIM can output 3*273119 matrices containing predicted signals for three histone modifications (H3K4me1, H3K4me3, and H3K36me3) at 273119 genomic regions when provided with a 3*273119 matrix containing H3K27ac, H3K27me3, and H3K9me3 signals. These matrices should be transformed into formatted read count matrices containing four columns (chr, start, end, signal).
(3) Predicted read count matrices were then converted to bedGraph format and utilized for subsequent analyses.
