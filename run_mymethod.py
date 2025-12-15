import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from sklearn.neighbors import LocalOutlierFactor

# -------------------------- Core Utility Functions --------------------------
# Get normalized count matrix (only keep outriderAE normalization as it's the core one used)
def get_cts(X, normalization):
    if normalization == 'outriderAE':
        cts = X['outriderNormCts']
    else:
        raise ValueError(f'Invalid normalization type: {normalization}')
    return cts

# Calculate z-score (only keep count-based z-score as core logic)
def zscore(df):
    df_zscore = df.transform(lambda x: (x - x.mean()) / x.std(), axis=1)
    return df_zscore.astype(np.float64)

# Calculate correlation matrix (only keep pearson correlation as core logic)
def corrmat(df):
    # Compute pearson correlation matrix
    corr = np.corrcoef(df)
    
    # Handle zero values to avoid calculation errors
    mask = (corr == 0)
    if mask.sum() > 0:
        corr[mask] = np.random.uniform(low=0, high=1e-3, size=mask.sum())
    
    corr = pd.DataFrame(corr, index=df.index, columns=df.index).astype(np.float64)
    np.fill_diagonal(corr.values, 0)  # Set diagonal to 0 (self-correlation)
    return corr

# Calculate devi index (core feature for AXO)
# abs_func: absolute value function
# corr: correlation matrix
# df_score: score matrix (outriderP)
# window_size: top N correlated genes
# g: current gene
def devi(abs_func, corr, df_score, window_size, g):
    g_corr = corr[g]
    z_self = df_score.loc[g, :]
    
    # Get top N positively correlated genes
    top_pos = g_corr.nlargest(window_size).index
    mean_z_pos = df_score.loc[top_pos, :].mean()
    
    # Only use positive correlation (mode 1, core logic of AXO)
    return abs_func(mean_z_pos - z_self)

# Outlier detection with LOF (core model for AXO)
def compute_prediction(X, n_neighbors=20):
    clf = LocalOutlierFactor(n_neighbors=n_neighbors, contamination="auto")
    clf.fit(X)
    return clf.negative_outlier_factor_

# -------------------------- Main Process --------------------------
def main():
    if len(sys.argv) != 4:
        print("Usage: python run_mymethod.py <cout_outsingle> <out_outrider> <out_mymethod>")
        sys.exit(1)
    
    # Get input parameters
    out_outsingle, out_outrider, out_mymethod = sys.argv[1], sys.argv[2], sys.argv[3] # {outsingle_out} {outrider_out} {mymethod_out}
    # Load core features (outsingleP and outrider results)
    X = {}
    # Load outsingle P-value (convert to -log)
    X['outsingleP'] = -np.log(pd.read_csv(out_outsingle, sep='\t', index_col=0))
    # Load outrider P-value (convert to -log)
    X['outriderP'] = -np.log(pd.read_csv(out_outrider, sep='\t', index_col=0))
    # Load outrider normalized counts and z-score (required for core calculation)
    X['outriderNormCts'] = pd.read_csv(f'{out_outrider}_normcounts.gz', sep='\t', index_col=0)
    X['outriderCtsZ'] = pd.read_csv(f'{out_outrider}_zScore.gz', sep='\t', index_col=0)
    
    # Core parameters (simplified, only keep the necessary ones for AXO)
    normalization_type = 'outriderAE'
    percent_choice = 0.1  # Top 10% correlated genes
    window_size = int(X['outriderNormCts'].shape[0] * percent_choice)
    score_type = 'outriderP'
    
    # Calculate devi feature (core step for AXO)
    cts = get_cts(X, normalization_type)
    df_zscore = zscore(cts)
    corr = corrmat(df_zscore)
    genes = df_zscore.index
    # Calculate devi for each gene
    res_list = [devi(abs, corr, X[score_type], window_size, g) for g in genes]
    X['devi_feature'] = pd.concat(res_list, axis=1, keys=genes).T
    
    # Feature combination: outsingleP + devi_feature (core combination for AXO)
    feature_combination = [X['outsingleP'], X['devi_feature']]
    Xsmall = pd.concat([f.stack() for f in feature_combination], axis=1, keys=['outsingleP', 'devi_feature'])
    
    # Predict with LOF (core model)
    array_g = X['outriderP'].index
    
    res_in_rows = [compute_prediction(Xsmall.loc[g, :].values) for g in array_g]
    df_res = pd.DataFrame(res_in_rows, index=array_g, columns=X['outriderP'].columns).astype(np.float64)
    
    # Save results
    outfile = f'{out_mymethod}' #_res_f0.1L20.txt'
    df_res.round(2).to_csv(outfile, sep='\t')
    os.system(f'gzip -f {outfile}')
    print(f"Result saved to: {outfile}.gz")

if __name__ == "__main__":
    main()