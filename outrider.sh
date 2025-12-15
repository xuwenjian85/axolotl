# outrider.sh

# 20250924 install by conda (https://anaconda.org/bioconda/bioconductor-outrider)
# conda install bioconda::bioconductor-outrider

incts=$1
norm_cts=${incts}.otrd_norm_cts.txt
out_otrd=${incts}.otrd.txt
out_otrd_gz=$2

demo_path=/media/eys/xwj/axolotl_rev/axo
data_path=`dirname ${incts}`

uid=`id | awk '{match($0,/([0-9]+)/,a)} {print a[0]}'`
echo ${uid}

# source /public/home/test1/soft/anaconda3/etc/profile.d/conda.sh; \
source /mnt/disk7t/xwj/soft/miniforge3/etc/profile.d/conda.sh; \
conda activate r4_outrider; \

Rscript --version; \
Rscript /mnt/disk7t/xwj/axolotl_rev/axo/script/outrider_pred.R \
${incts} \
${norm_cts} \
${out_otrd}

# Normalize raw counts using Outrider
# docker run -u $uid:$uid -v ${demo_path}:/axo -v ${data_path}:/data  \
# --rm -it r4.2:jammy bash -c \
# "Rscript /axo/script/outrider_pred.R \
# /data/`basename ${incts}` \
# /data/`basename ${norm_cts}` \
# /data/`basename ${out_otrd}`"

gzip -c ${out_otrd} > ${out_otrd_gz}

# docker save -o r4.2_jammy.tar r4.2:jammy
# chmod 755 r4.2_jammy.tar