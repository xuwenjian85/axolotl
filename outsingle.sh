# outsingle.sh

incts=$1
out_osg=${incts}.osg.txt
out_osg_gz=$2

uid=`id | awk '{match($0,/([0-9]+)/,a)} {print a[0]}'`
echo "$uid"

demo_path=/media/eys/xwj/axolotl_rev/axo
data_path=`dirname ${incts}`

# source /public/home/test1/soft/anaconda3/etc/profile.d/conda.sh; \
# conda activate OutSingle

python /mnt/disk7t/xwj/axolotl_rev/axo/script/outsingle.py \
${incts} \
${out_osg} \
/mnt/disk7t/xwj/axolotl_rev/axo/script/outsingle/

# 2. Calculate aberrant p-values using OutSingle
# docker run -u $uid:$uid -v ${demo_path}:/axo -v ${data_path}:/data   \
# --rm -it py3 bash -c \
# "python /axo/script/outsingle.py \
# /data/`basename ${incts}` \
# /data/`basename ${out_osg}` \
# /axo/script/outsingle/"

gzip -c ${out_osg} > ${out_osg_gz}

# docker save -o py3.tar py3
# chmod 755 py3.tar
