# outsingle_helper.py
# input: 
#   summary table & job arguments
# output: 
#   file paths for outsingle script

import sys
import os
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

import pandas as pd
import numpy as np
  
d_sum = pd.read_csv(sys.argv[1], sep="\t", index_col=0, header=0)
dataset = d_sum.index[ np.int32(sys.argv[2]) ]
print( sys.argv, dataset)
# file paths 
# cts_col = 'NewIntegerCounts' if 'NewIntegerCounts' in d_sum.columns else 'RawCounts'
cts_col = 'RawCounts'
outscore_col ='outsingleRes'
# read job from df_summary & job_id
file_incts = d_sum.loc[dataset, cts_col] 
file_outscore = d_sum.loc[dataset, outscore_col]

OutSingle='/media/eys/xwj/expression_outlier/run_outsingle.py'
print( OutSingle, file_incts, file_outscore)
os.system('python %s %s %s' % (OutSingle, file_incts, file_outscore))