# Ensembles of genome-coverage single-cell histone modifications reveal epigenetic lineages during mouse preimplantation development

## Brief Descriptions of Analysis Scripts

### Processing TACIT data
Use tacit.sh in the `TACIT` folder to analyze raw sequencing data, including mapping to the mm10 reference genome and removing PCR duplicates.

### Processing CoTACIT data
1. Split raw sequencing data for each histone modification based on T5 and T7 barcodes with scripts in the `CoTACIT` folder.
2. Analyze the split sequencing data for each histone modification using tacit.sh.

### Processing scRNA-seq data
1. Split raw sequening data for each single cell based on oligo dT barcodes with scripts in the `scRNA-seq/00_split_single_cells` folder.
2. Analyze sequencing data including cutting adaptor, mapping to mm10 reference genome, and calculating TPM, with scripts in the `scRNA-seq/01_qc` folder.

### Integration analyses
1. Multimodality integration with scRNA-seq as anchors follows the Seurat V4 anchoring workflow: https://satijalab.org/seurat/. Detailed scripts are located in the `fig2` folder.
2. Pipeline for WNN analysis: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis. Detailed scripts are located in the `fig2` and `fig3` folder.

### ChromHMM analyses
1. ChromHMM manual: https://ernstlab.biolchem.ucla.edu/ChromHMM/
2. Filter out all reads from mitochondrial DNA to reduce noise and interference.
3. For Fig.4, all ChromHMM model and emission probabilities txt file for synthetic cells are available in the `fig4/data` folder. And genome-wide chromatin state annotations at a 200-kb resolution for all cells are provided. Detailed scripts are located in the `fig4` folder.
4. For Fig.5, the forward-backward algorithm is used to learn the posterior probability distribution for all synthetic cells, with the model text file provided in the `fig5/data` folder. Genome-wide chromatin state annotations at a 2000-kb resolution for all cells are also included. Detailed scripts are located in the `fig5` folder.

## Step-by-step codes for reproducing each main figure are available in their respective subfolders.

