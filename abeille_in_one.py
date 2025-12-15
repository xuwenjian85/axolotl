# abeille in one script 20251111

# usage: 
#   conda activate py3.10_tf2.8
#   python abeille_in_one.py 

import sys
import os
import pandas as pd
import subprocess
import yaml

# 使用时改为本地路径
ctsFile = sys.argv[1] # input文件路径，基因表达量
out_abeille_gz = sys.argv[2] # output文件路径，打分矩阵

# 读取配置文件20251111
current_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(current_dir, 'config.yaml'), 'r') as f:
    config = yaml.safe_load(f)
conda_path = config['paths']['conda']
# conda_path = "/mnt/disk7t/xwj/soft/miniforge3/etc/profile.d/conda.sh" # 本地conda.sh所在目录 
script_path = config['paths']['script1']
# script_path = "/mnt/disk7t/xwj/axolotl_rev/script" # abeille py和r脚本所在文件夹
r_conda_env = config['paths']['abeille_env_r'] #'r4.0_abeille' # 如有必要，改为你R环境的名称 
# abeille_env_r: r4.0_abeille
# ---------
# 以下内容保持不变
script = os.path.join(script_path, 'abeille_cpu_pred.R')
functionfile = os.path.join(script_path, 'IdentifyAGE_no_selection.R')
# 添加 subfolder 到 sys.path
sys.path.append(os.path.join(script_path))

############VAE in 当前python环境
from abeille_tab import abeille_VAE

ReconsCtsFile = f'{ctsFile}.ReconsCts'
ResFile = f'{ctsFile}.abeille'

print(ctsFile, ReconsCtsFile, ResFile)

cts_recons = abeille_VAE(ctsFile)
cts_recons.to_csv(ReconsCtsFile, sep="\t", index=True, header=True)

###########IdentifyAGE in R

conda_cmds = [
    f'source {conda_path}',
    f'conda activate {r_conda_env}',
    'Rscript --version',
    f'Rscript {script} {ctsFile} {ReconsCtsFile} {ResFile} {functionfile}',  
    f'gzip -c {ResFile} > {out_abeille_gz}', 
]
# 运行
subprocess.run(['bash', '-c', '; '.join(conda_cmds)])
