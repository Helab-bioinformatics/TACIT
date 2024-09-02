import pandas as pd
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr

def calculate_metrics(sample_file, true_file):
    pred_df = pd.read_csv(sample_file, sep='\s+', header=None)
    true_df = pd.read_csv(true_file, sep='\s+', header=None)
    
    # Assume the columns correspond to H3K4me3, H3K4me1, H3K36me3 respectively
    metrics = {}
    for i, histone_mark in enumerate(['H3K4me3', 'H3K4me1', 'H3K36me3']):
        y_pred = pred_df[i]
        y_true = true_df[i]
        
        y_pred = y_pred.astype(float)
        y_true = y_true.astype(float)
        
        # Compute Mean Squared Error (MSE)
        mse = mean_squared_error(y_true, y_pred)
        
        # Compute Pearson correlation
        corr, _ = pearsonr(y_true, y_pred)
        
        metrics[histone_mark] = {
            'MSE': mse,
            'Correlation': corr
        }
    
    return metrics

sample_file = 'path1'  #Path for prediction files (files needs to be transposed)
true_file = 'path2' #Path for ground truth files (files needs to be transposed)

metrics = calculate_metrics(sample_file, true_file)

for histone_mark, metric in metrics.items():
    print(f"{histone_mark}: MSE = {metric['MSE']}, Correlation = {metric['Correlation']}")
