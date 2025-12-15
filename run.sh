
task=0
task_config=$1
script=parse_task_abeille.py
source /mnt/disk7t/xwj/soft/miniforge3/etc/profile.d/conda.sh
conda activate py3.12
# python ${script} ${task_config} ${task}

###
task=0
task_config=$1
script=parse_task_outrider.py
source /mnt/disk7t/xwj/soft/miniforge3/etc/profile.d/conda.sh
conda activate py3.12
python ${script} ${task_config} ${task}

###
task=0
task_config=$1
script=parse_task_outsingle.py
source /mnt/disk7t/xwj/soft/miniforge3/etc/profile.d/conda.sh
conda activate py3.12
python ${script} ${task_config} ${task}

### 
task=0
task_config=$1
script=parse_task_axo.py
source /mnt/disk7t/xwj/soft/miniforge3/etc/profile.d/conda.sh
conda activate py3.12
python ${script} ${task_config} ${task}
