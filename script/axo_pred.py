## running time is depending on the Sample numbers
# 15min(Ns=100), 20 min(Ns=150) , 40 min(Ns=300), 50min(Ns=423) ~ 1 hr,

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
import time, sys
from multiprocessing import Pool
from copy import deepcopy
import pandas as pd
import numpy as np
from sklearn import metrics
from sklearn.covariance import EllipticEnvelope
from sklearn.svm import OneClassSVM
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor

### helper functions
## correlations of gene pair
def pearsons(data):
    pcc = pd.DataFrame(data=metrics.pairwise.cosine_similarity(data, dense_output=True), index = data.index, columns = data.index)
    if (pcc == 0).sum().sum() > 0:
        zero_x, zero_y = np.where(pcc == 0)
        for r, c in zip(zero_x, zero_y):
            pcc.iloc[r,c] = np.random.uniform(low=0, high=1e-3)
        assert (pcc == 0).sum().sum() == 0
    np.fill_diagonal(pcc.values, 0)
    return pcc

def worker(sample, x,  zscore):     # usage: worker(  sample=df_cts.columns[0], x='z_rank', zscore =z_rank )
    tmp = pd.concat( [ zscore[ sample ].groupby( bin_array_dict[x][g] ).mean() for g in df_cts.index ], axis=1 )
    tmp = tmp.T
    tmp.index = df_cts.index
    col = -(tmp.subtract(zscore[sample], axis=0) )
    # don't **2
    return col

def devi(sample): 
    ## ONE sample,根据2种pcc计算devi & save results                
    # print(time.ctime(),'qbin=%s : %s' % (q, sample))
    devi_rank_df = worker( sample= sample, x='z_rank', zscore = z_rank )
    devi_logcts_df = worker( sample= sample, x='z_cts', zscore = z_cts )
    return devi_rank_df, devi_logcts_df

## 
print('\n'.join(sys.argv))
q = 50
df_cts = pd.read_csv(sys.argv[1], sep="\t", skiprows=0, index_col=0)
Xy = sys.argv[2]
ResMyModels=sys.argv[3]
ods_res = pd.read_csv(sys.argv[4], sep="\t", skiprows=0, index_col=0)
osg_res = pd.read_csv(sys.argv[5], sep="\t", index_col=0, header=0)
nCPU = np.int32(sys.argv[6])

HC = pd.DataFrame({
    'rank_mean':  df_cts.rank(axis='index').mean(axis='columns'), 
    'rank_std': df_cts.rank(axis='index').std(axis='columns') + 1e-3,
})

# Z_d_rank
z_rank = df_cts.rank().subtract( HC['rank_mean'], axis='index').divide(HC['rank_std'], axis='index')
z_rank = np.clip (z_rank , -20, 20).astype(np.float16)

## high abund genes should reduce z_rank
if True: # patch20221210 
    scale = [ np.int16(g+ 1)/100  if g < 100 else 1 for g in HC['rank_std'].values]
    z_rank = z_rank.multiply( scale, axis = 'index')

## Z_d_foldchange
z_cts = df_cts.transform(lambda x: (x-x.mean())/x.std(), axis=1)
z_cts = np.clip (z_cts , -20, 20).astype(np.float16)

print(z_cts.iloc[:3, :3] ) # index is ok, not changed

## two methods to pearson's Z(z_rank), Z(log, foldchange)
pcc = dict.fromkeys(['z_rank', 'z_cts'])
pcc['z_rank'] = pearsons(z_rank).astype(np.float16)
pcc['z_cts'] = pearsons(z_cts).astype(np.float16)    

# gene as column，other gene split to q bins by correlations values
bin_array_dict = dict.fromkeys(pcc) 
for x in bin_array_dict:
    bin_array_dict[x] = pd.concat( [ pd.qcut(pcc[x][g], q,  labels=False,  retbins=False) for g in pcc[x].index ], axis = 1 ).astype(np.int16)
    
print(time.ctime(), 'pcc done')
# deviation features for genes x samples
print(time.ctime(),'devi loop start')
res_in_rows = []
with Pool(nCPU) as pool:
    print(time.ctime(),'Pool map CPU=%s qbin=%s' % (nCPU, q))
    res_in_rows = pool.map( devi, df_cts.columns)

    
tmp=[]
for i in range(df_cts.shape[1]):
    s = df_cts.columns[i]
    res_scaled = deepcopy(res_in_rows[i])
    feature_one = pd.concat( res_scaled , axis =1, ignore_index=True )
    feature_one.index = pd.MultiIndex.from_product([feature_one.index, [s] ], names=['Gene','Sample'])
    tmp.extend([feature_one])
print(time.ctime(), 'devi done')

# with open(d_sum.loc[dataset, 'Xy']+'tmp', 'wb') as f:   
#     pickle.dump( tmp, f )

X = pd.concat(tmp)
X = X.loc[:, [0,50]].rename(columns={0:'rank_devi',50:'cts_devi'})
print(X.head())

# X.loc[:, 'rank_z'] = z_rank.loc[:, df_cts.columns ].stack()
# X.loc[:, 'cts_z'] = z_cts.loc[:, df_cts.columns ].stack()

y = pd.Series(data = 0, index = X.index, name='label')
# for idx, row in df_corrupt.iterrows():
#     y.loc[(row.corruptGene, row.simu )] = 1

# 20230601 修改保证y和X可以match, 测试
# for idx,gs in df_corrupt.iterrows():
#     eout_gs = (gs['corruptGene'],gs['simu'])
#     if eout_gs in y.index:
#         y.loc[eout_gs]  =1 

X = pd.concat([X,
               z_rank.stack().rename('rank_z'),
               z_cts.stack().rename('cts_z'),
               y], axis=1)

X.to_csv(Xy)

##### 
# X = pd.read_pickle(d_sum.loc[dataset, 'Xy'])
X = X.dropna()

## take -log(pvalues) as my input features
X.loc[:, ['ods_pv' ] ] = (ods_res.rename(columns={'geneID':'Gene', 'sampleID': 'Sample'})
    .set_index(['Gene','Sample']).loc[:, ['pValue']].copy()
    .loc[X.index].values)
X.loc[:, ['ods_pv' ] ] = - X.loc[:, ['ods_pv'] ].transform(np.log).values

X.loc[:, ['osg_pv' ] ] = (osg_res.melt(ignore_index=False, var_name = ['Sample'],value_name='ds')
        .set_index(['Sample'],append=True).loc[X.index].values)
X.loc[:, ['osg_pv' ] ] = - X.loc[:, ['osg_pv'] ].transform(np.log).values

dict_model = {
    'L': 'LOF', 
    # 'I': 'iForest',
    # 'O': 'OCSVM',
}

dict_feature = {
    'f-0' : ['cts_z','rank_z','rank_devi','cts_devi','osg_pv'],
    # 'f-12': [                 'rank_devi','cts_devi','osg_pv'],
    # 'f-34': ['cts_z','rank_z',                       'osg_pv'],
    # 'f-5' : ['cts_z','rank_z','rank_devi','cts_devi',        ],
}

## save result
dict_res = dict.fromkeys([ (f+m) for f in dict_feature for m in dict_model ]) 

# Outlier detection: single gene per loop
def compute_prediction(X, model_name):

    # print(f"Computing {model_name} prediction...")
    if model_name == "EEnv": # One-class classification # 6min 38s
        clf = EllipticEnvelope(contamination=0.001, support_fraction = 0.99 )
        y_pred = clf.fit(X).decision_function(X)
        
    if model_name == "LOF": #"distance-based"# <1min
        clf = LocalOutlierFactor(n_neighbors=20, contamination="auto")
        clf.fit(X)
        y_pred = clf.negative_outlier_factor_
        
    if model_name == "iForest": #"distance-based" # 1min
        clf = IsolationForest(n_estimators=20, contamination='auto', max_samples =1.0) # 100,20;contamination=0.01
        y_pred = clf.fit(X).decision_function(X)
        
    if model_name == 'OCSVM': # One-class classification
        clf = OneClassSVM(gamma='auto')
        y_pred = clf.fit(X).decision_function(X)

    return y_pred

def parallel_prediction(g):
    y_pred = compute_prediction(Xsmall.loc[g,:].values, model_name)
    return y_pred

# one df result from one combination setting of (feature, model )
array_g = X.index.get_level_values(level = 0).unique()
array_s = X.index.get_level_values(level = 1).unique()

for feature_code in dict_feature:
    for model_code in dict_model:
        print(feature_code, dict_feature[feature_code], model_code, dict_model[model_code])
        Xsmall = X[ dict_feature[ feature_code ] ].copy() 
        model_name = dict_model[ model_code ] 
        
        # 1. model fit
        with Pool(nCPU) as pool: 
            print(time.ctime(), 'Pool map ->',nCPU)
 
            res_in_rows = pool.map( parallel_prediction, array_g)
        # my result
        dict_res[ feature_code + model_code ] = pd.DataFrame(res_in_rows, index=array_g, columns=array_s )
    
# save
f,m = ['f-0','L']
    
dict_res[f+m].to_csv(ResMyModels)
