date
hostname

nCPU=5

incts=testdata/df_cts.txt

out_osg=${incts}.osg.txt
norm_cts=${incts}.otrd_norm_cts.txt
out_otrd=${incts}.otrd.txt
Xy=${incts}.axo_xy.txt
out_axo=${incts}.axo.txt

OutSingleBin='/home/xwj/work/axo/script/outsingle/'
uid=`id | awk '{match($0,/([0-9]+)/,a)} {print a[0]}'`
echo "$uid"

# 1. Normalize raw counts using Outrider
docker run -u $uid:$uid -v `pwd`:/axo \
--rm -it r4.2:jammy_outrider bash -c \
"Rscript /axo/script/outrider_pred.R /axo/${incts} /axo/${norm_cts} /axo/${out_otrd}"

# 2. Calculate aberrant p-values using OutSingle
docker run -u $uid:$uid -v `pwd`:/axo \
--rm -it py3 bash -c \
"python /axo/script/outsingle.py /axo/${incts} /axo/${out_osg} /axo/script/outsingle/;\
wc /axo/script/outsingle.py /axo/${incts} /axo/${out_osg}"

# 3. Pred aberrant score with axolotl
docker run -u $uid:$uid -v `pwd`:/axo \
--rm -it py3 bash -c \
"python /axo/script/axo_pred.py \
/axo/${norm_cts} /axo/${Xy} /axo/${out_axo} /axo/${out_otrd} /axo/${out_osg} $nCPU"