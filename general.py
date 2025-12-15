import pandas as pd
import numpy as np
from sklearn.metrics import auc

# 定义一个函数，将大于0.01的数值四舍五入到三位小数
def round_three_decimal(x):
    if x > 0.01:
        return round(x, 3)
    return x

# subset from full cts, 返回data
def make_cts( cts, samples, simu):
    cols = samples.loc[simu,:]
    cts_subset = cts.loc[:, cols ]
    return cts_subset

# 各种方法原始输出里包含中间结果、长浮点数，保存更少的有效数字，从而节省空间
def postprocess(method,file, file_cts):
    outfile = file
    print(method,file)
    cts = pd.read_csv(file_cts, sep='\t',index_col='Gene')

    if method == 'outrider':
        ################# outrider
        # file = sys.argv[1]
        outfile = file
        df = pd.read_csv(file, sep='\t',index_col=0)
        c = 'raw' # 先保留一个原结果
        df.to_csv(f'{outfile}_{c}.gz', sep='\t')
        df = df.rename(columns={'geneID':'Gene', 'sampleID': 'Sample', 'padjust':'ds'})
        # --- other values
        # zScore
        # l2fc
        # normcounts
        for c in ['zScore', 'l2fc','normcounts' ]:
            # df.pivot(index = 'Gene', columns='Sample',values = c) \
            #     .loc[cts.index, cts.columns]\
            #     .to_csv(f'{outfile}_{c}.gz', sep='\t')    
            df_wide = df.pivot(index = 'Gene', columns='Sample',values = c) 
            common_genes = [gene for gene in cts.index if gene in df_wide.index]
            df_wide.loc[common_genes, cts.columns].to_csv(f'{outfile}_{c}.gz', sep='\t')
        # --- pvalue
        c = 'pValue'
        # df.pivot(index = 'Gene', columns='Sample',values = c ) \
        #     .applymap(round_three_decimal)\
        #     .loc[cts.index, cts.columns]\
        #     .to_csv(f'{outfile}_{c}.gz', sep='\t')
        df_wide = df.pivot(index = 'Gene', columns='Sample',values = c) \
                    .applymap(round_three_decimal)
        common_genes = [gene for gene in cts.index if gene in df_wide.index]
        df_wide.loc[common_genes, cts.columns].to_csv(f'{outfile}_{c}.gz', sep='\t')
                   
        # --- padjust
        # df.pivot(index = 'Gene', columns='Sample',values = 'ds') \
        #     .applymap(round_three_decimal)\
        #     .loc[cts.index, cts.columns]\
        #     .to_csv(outfile, sep='\t')
        df_wide = df.pivot(index = 'Gene', columns='Sample',values = 'ds') \
                    .applymap(round_three_decimal)
        common_genes = [gene for gene in cts.index if gene in df_wide.index]
        df_wide.loc[common_genes, cts.columns].to_csv(outfile, sep='\t')
            
    elif method == 'outsingle':   
        # outsingle_score = outsingle_score.applymap(round_three_decimal)
        df = pd.read_csv(file, sep='\t',index_col=0)
        # 使用applymap函数应用该函数到DataFrame的每个元素
        # df = df.applymap(round_three_decimal).loc[cts.index, cts.columns]
        # df.to_csv(outfile, sep='\t')
        df = df.applymap(round_three_decimal)
        common_genes = [gene for gene in cts.index if gene in df.index]
        df.loc[common_genes, cts.columns].to_csv(outfile, sep='\t')
        
        
    elif method == 'abeille' :
        # abeille
        c = 'raw' # 先保留一个原结果
        df = pd.read_csv(file, sep='\t',index_col=0)
        df.to_csv(f'{outfile}_{c}.gz', sep='\t')
        # 再后续处理
        # df = df.rename(columns={'Transcript': 'Gene', 'anomaly_score':'ds'})\
        #     .pivot(index = 'Gene', columns='Sample',values = 'ds') \
        #     .applymap(round_three_decimal)\
        #     .loc[cts.index, cts.columns]
        # df.to_csv(outfile, sep='\t')
    else:
        pass
    
# 函数定义

## precision recall curve add segments for easy plotting
# n_seg = 100
# recall_seg = np.linspace(0, 1,num=n_seg+1)
def pred_outlier_auprc(df_pred, df_outlier, n_seg=100):
    
    df_pred.index.name = 'Gene'
    df_pred.columns.name = 'Sample'
    pr_curve = df_pred.T.melt(ignore_index=False,value_name='pred').set_index(['Gene'],append=True).reorder_levels([1, 0])
    pr_curve = pr_curve.sort_values(by='pred',ascending=True)
    pr_curve['label'] = 0
    
    for i, row in df_outlier.iterrows(): 
    #     print(i,row['Gene'],row['Sample'])
        pr_curve.loc[(row['Gene'], row['Sample']), 'label' ] = 1
    pr_curve['rank'] = range(1,pr_curve.shape[0]+1,1)
    pr_curve['precision'] = (pr_curve['label'].cumsum()/pr_curve['rank']).round(4)
    pr_curve['recall'] = (pr_curve['label'].cumsum()/pr_curve['label'].sum()).round(4)

    # auprc value
    auprc = auc(pr_curve.loc[:,'recall'],pr_curve.loc[:,'precision'])
    
    # save a simple curve for plotting
    pr_curve = pr_curve.drop_duplicates(subset=['precision', 'recall'])
    # keep equally distributed points
    recall_seg = np.linspace(0, 1,num=n_seg+1)
    idx = [ (pr_curve['recall'] - seg).abs().idxmin() for seg in recall_seg ] # type: ignore
    
    return auprc, pr_curve.loc[idx,]