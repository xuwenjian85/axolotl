source /mnt/disk7t/xwj/soft/miniforge3/etc/profile.d/conda.sh
conda activate py3.10_tf2.8
ctsfile=$1
outscore=$2
python /mnt/disk7t/xwj/axolotl_rev/script/abeille_in_one.py ${ctsfile} ${outscore}