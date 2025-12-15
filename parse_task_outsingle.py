import sys
import os
import pandas as pd
import numpy as np
import tempfile
import os
from general import make_cts, postprocess

# 接收job
task_config_path=sys.argv[1] 
jobid=np.int32(sys.argv[2])

# 任务列表
task_config = pd.read_csv(task_config_path, sep='\t')

# GTEx one tissue全部样本的表达量矩阵
cts = pd.read_csv(task_config.loc[jobid,'cts'], sep='\t',index_col='Gene')
# 指定seed随机抽样得到的样本子集
samples = pd.read_csv(task_config.loc[jobid,'samples'], sep='\t',index_col='task')
print(cts.shape, samples.shape, jobid)
# 任务编号，对应一个样本子集
simu = jobid

# 临时目录，存储各种中间文件, 结束即被自动删除
with tempfile.TemporaryDirectory() as tmp_dir:
    print('临时目录:', tmp_dir)
    cts_tmpfile = f'{tmp_dir}/{simu:03d}.txt'
    cts_subset = make_cts(cts, samples, simu)
    cts_subset.to_csv( cts_tmpfile, sep='\t')

    # OUTRIDER
    
    # OUTSINGLE
    script_outsingle = '/mnt/disk7t/xwj/axolotl_rev/script/outsingle.sh'
    outscore = task_config.loc[simu,'OUTSINGLE']
    os.system(f'bash {script_outsingle} {cts_tmpfile} {outscore}')
    postprocess('outsingle', outscore, cts_tmpfile)
    
    # ABEILLE
    
    # MY
    os.system(f'tree -sh {tmp_dir}') # checkpoint
