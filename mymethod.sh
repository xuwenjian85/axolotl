# mymethod.sh
# call run_mymethod.py 
# 
source /mnt/disk7t/xwj/soft/miniforge3/etc/profile.d/conda.sh
conda activate py3.12

script=/mnt/disk7t/xwj/axolotl_rev/script/run_mymethod.py
outsingle_output=$1
outrider_output=$2
mymethod_output=$3
python ${script} ${outsingle_output} ${outrider_output} ${mymethod_output}