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
1. Multimodality integration with scRNA-seq as anchors follows the Seurat V4 anchoring workflow: https://satijalab.org/seurat/.
2. WNN analysis: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.
2. scDeepHIM:
