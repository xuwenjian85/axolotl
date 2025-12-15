
# gtex 49 tissues
for file in '/mnt/disk7t/xwj/axolotl_rev/result/dataset_gtex_nmd/task_config/*.config'; do
    bash run_abeille.sh "$file"
    bash run_outsingle.sh "$file" 
    bash run_outrider.sh  "$file"
    bash run_axo.sh  "$file"
done
#-------------------------------- --------------------------------------------------------------------------
# kremer119 14k all samples
for file in '/mnt/disk7t/xwj/axolotl_rev/result/dataset_gtex_nmd/task_config/t00_FB_s119_g14409.config'; do
    bash run_abeille.sh "$file"
    bash run_outsingle.sh "$file" 
    bash run_outrider.sh  "$file"
    bash run_axo.sh  "$file"
done

#-------------------------------- --------------------------------------------------------------------------
# kremer119 11k size 30~100 with different seeds
for file in /mnt/disk7t/xwj/axolotl_rev/result/dataset_nc2017kremer/task_config/t01_FB_s119_g11048_nmd_size*.config; do

    bash run_abeille.sh "$file"
    bash run_outsingle.sh "$file" 
    bash run_outrider.sh  "$file"
    bash run_axo.sh  "$file"
done

for file in /mnt/disk7t/xwj/axolotl_rev/result/dataset_nc2017kremer/task_config/t00_FB_s119_g11048_pct*.config; do
    bash run_abeille.sh "$file"
    bash run_outsingle.sh "$file" 
    bash run_outrider.sh  "$file"
    bash run_axo.sh  "$file"
done


#-------------------------------- --------------------------------------------------------------------------
# pmuscle36 13.5k genes all samples
for file in /mnt/disk7t/xwj/axolotl_rev/result/dataset_pmuscle_36/task_config/t00_M_s36_g13573.config; do
    bash run_abeille.sh "$file"
    bash run_outsingle.sh "$file" 
    bash run_outrider.sh  "$file"
    bash run_axo.sh  "$file"
done
# pmuscle36 12k genes all samples with different seeds
for file in /mnt/disk7t/xwj/axolotl_rev/result/dataset_pmuscle_36/task_config/t00_M_s36_g13573_size12000_seed*.config; do
    bash run_abeille.sh "$file"
    bash run_outsingle.sh "$file" 
    bash run_outrider.sh  "$file"
    bash run_axo.sh  "$file"
done

#-------------------------------- --------------------------------------------------------------------------
# pfib ss 269 samples 
for file in /mnt/disk7t/xwj/axolotl_rev/result/dataset_pfib_423_split/task_config/t00_FBSS_s269_g12369.config; do
    bash run_abeille.sh "$file"
    bash run_outsingle.sh "$file" 
    bash run_outrider.sh  "$file"
    bash run_axo.sh  "$file"
done

# pfib ns 154 samples 
for file in /mnt/disk7t/xwj/axolotl_rev/result/dataset_pfib_423_split/task_config/t01_FBNS_s154_g13411.config; do
    bash run_abeille.sh "$file"
    bash run_outsingle.sh "$file" 
    bash run_outrider.sh  "$file"
    bash run_axo.sh  "$file"
done

# size 50 or 100 at different outlier pct
for file in /mnt/disk7t/xwj/axolotl_rev/result/dataset_pfib_423_split/task_config/t0*_FB*pct*size*.config; do
    bash run_abeille.sh "$file"
    bash run_outsingle.sh "$file" 
    bash run_outrider.sh  "$file"
    bash run_axo.sh  "$file"
done


#-------------------------------- --------------------------------------------------------------------------

for file in /mnt/disk7t/xwj/axolotl_rev//result/dataset_ly111/task_config/t00_B_s111_g11310.config; do
    bash run_abeille.sh "$file"
    bash run_outsingle.sh "$file" 
    bash run_outrider.sh  "$file"
    bash run_axo.sh  "$file"
done