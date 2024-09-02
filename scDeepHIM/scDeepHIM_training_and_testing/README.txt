scDeepHIM Training and Testing Guide
Follow these step-by-step instructions to train and test scDeepHIM:

1. Generate Count Matrices:
Create read count matrices for histone modifications for each sample in the training group using the bamCoverage function (00_Data_pre-Treatment.sh) or an alternative method.
2. Prepare Input and Target Files:
Partition the matrices into input files (H3K27ac, H3K27me3, H3K9me3) and target files (H3K4me1, H3K4me3, H3K36me3).
Train the scDeepHIM models using profiles of six histone modifications, allocating 80% of the samples for training and 20% for evaluation.
Place the input (H3K27ac, H3K27me3, H3K9me3) and target (H3K4me1, H3K4me3, H3K36me3) files for both training and testing samples in separate directories.
Load these files during model training (01_scDeepHIM_Model_training_1000echo_GPU.py).
3. Load and Predict:
Load the files to be predicted into the model using the input files with H3K27ac, H3K27me3, and H3K9me3 signals (02_Prediction_scDeepHIM.py).
During training, the scDeepHIM model optimizes an objective function to predict H3K4me1, H3K4me3, and H3K36me3 signals from both bulk and single-cell ChIP-seq profiles.
After training, scDeepHIM outputs matrices containing predicted signals for the three histone modifications (H3K4me1, H3K4me3, H3K36me3) across 3,273,119 genomic regions.
4. Evaluate Prediction Performance:
Evaluate the model’s prediction performance using 03_Performance_evaluation.py.
5. Format Output Matrices:
Convert the predicted matrices into formatted read count matrices with four columns (chr, start, end, signal) using 04_Prediction_to_readcount.R.
6. Convert to bedGraph Format:
Transform the predicted read count matrices to bedGraph format for further analysis (05_Prediction_txt_to_bedgragh.py and 06_Prediction_bedgragh_to_bigwig.sh).
