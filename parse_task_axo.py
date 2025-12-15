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

    # ABEILLE
    
    # MY
    # 获取样本总数
    # 假设用dim_q0 值，已经运行完成OUTRIDER+dim_q0打分
    Nsample = cts_subset.shape[1]
    div_q0 = 3
    dim_q0 = int(Nsample / div_q0)
    print("q0_value:", div_q0)
    print(f'Running OUTRIDER with dim_q0 = {dim_q0}')
    # OUTRIDER default输出路径衍生出dim_q0输出路径
    outscore_OUTRIDER = task_config.loc[simu, 'OUTRIDER']
    outscore_with_dim = f'{outscore_OUTRIDER}_dim{dim_q0}.gz'
    
    if os.path.isfile(outscore_with_dim): # 如果文件存在则跳过计算20251112
        print(f'{outscore_with_dim}, file exists')
    else:
        print(f'{outscore_with_dim}, compute start.')
        # 读取配置文件
        import yaml
        import subprocess
        
        current_dir = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(current_dir, 'config.yaml'), 'r') as f:
            config = yaml.safe_load(f)
        # 使用路径
        conda_path = config['paths']['conda']
        script_path = config['paths']['script1']
        # 脚本路径
        script = os.path.join(script_path, 'outrider_pred_dim.R')
        # 准备命令, 无需脚本, sh直接调用
        conda_cmds = [
            f'source {conda_path}',
            'conda activate r4_outrider',
            'Rscript --version',
            f'Rscript {script} {cts_tmpfile} {dim_q0} {cts_tmpfile}.otrd_norm_cts.txt {cts_tmpfile}.otrd.txt',  
            f'gzip -c {cts_tmpfile}.otrd.txt > {outscore_with_dim}', 
        ]
        # 运行
        subprocess.run(['bash', '-c', '; '.join(conda_cmds)])
        # 检查点：打印临时目录的内容
        os.system(f'tree -sh {tmp_dir}')
        # 后处理
        postprocess('outrider', outscore_with_dim, cts_tmpfile)
        
    script = '/mnt/disk7t/xwj/axolotl_rev/script/mymethod.sh'
    outscore = task_config.loc[simu,'MyMethod'] # filename not shown q0=3
    outscore1 = task_config.loc[simu,'OUTSINGLE']
    outscore2 = outscore_with_dim
    
    print(outscore, outscore1, outscore2)
    os.system(f'bash {script} {outscore1} {outscore2} {outscore}')

    os.system(f'tree -sh {tmp_dir}') # checkpoint
