# helper script for OutSingle:
# MyHelperScript=/media/eys/xwj/expression_outlier/helper_outsingle.py

import sys
import os
import numpy as np
import pandas as pd

# read job from df_summary & job_id
incts=sys.argv[1] #'/media/eys/xwj/soft/outsingle/testdata/df_cts_muscle36.txt'
outscore=sys.argv[2]#'/media/eys/xwj/soft/outsingle/testdata/df_cts_muscle36.score.tsv'

# path to OutSingle
OutSingleBin=sys.argv[3]#'/media/eys/xwj/soft/outsingle/'

# Process input matrix, remove row names, and retain column names
# os.system('cut -f1- %s > %s' % (incts, incts+'_norowname.txt'))
temp = pd.read_csv(incts,sep="\t",index_col=0, header=0)
temp.index.name=None
temp.to_csv(incts + '_norowname.txt', sep='\t')

# z-score estimation; calculate the final OutSingle score
os.system('python %s/fast_zscore_estimation.py %s' % (OutSingleBin, incts+'_norowname.txt'))
os.system('python %s/optht_svd_zs.py %s' % (OutSingleBin, incts+'_norowname-fzse-zs.csv'))

# modify the OutSingle score table by adding the rownames of count table
rowname = pd.read_csv(incts, sep="\t", index_col=0, header=0)
outsingle_score = pd.read_csv(incts+'_norowname-fzse-zs-svd-optht-zs-pv.csv', 
                              sep="\t", index_col=0, header=0)

outsingle_score.index = rowname.index
print(rowname.shape, outsingle_score.shape)
outsingle_score.to_csv(outscore, sep="\t", index=True, header=True)

os.system('rm %s' % incts+'_norowname*')