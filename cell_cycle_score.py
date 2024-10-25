
import os,sys
import scanpy as sc
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
import matplotlib.pyplot as plt
import celltypist
from celltypist import models
sc.settings.set_figure_params(dpi=300,dpi_save=300, facecolor='white')

adata1 = sc.read_h5ad('~/scABE_filter_anno_day31421_fullgene_public.h5ad')
sc.pp.normalize_total(adata1, target_sum=1e4)
sc.pp.log1p(adata1)
adata1

genelist = pd.read_excel('~/genelist.xlsx')#,sheet_name='Gene Sets Used in Analysis'
counts = adata1.to_df().T
df = pd.merge(genelist,counts,left_on='gene_name',right_on=counts.index,how='inner')


def genelist_filter(phase):
    genelist_S = df[df['phase']==phase]
    counts_S = pd.merge(genelist_S[['phase','gene_name']],counts,left_on='gene_name',right_on=counts.index,how='left')
    counts_S.index = counts_S['gene_name']
    counts_S = counts_S.iloc[:,2:]
    counts_S.loc['mean',:] = counts_S.mean()
    counts_S = counts_S.loc['mean',]
    corr_S = pd.DataFrame()
    corr_S['S'] = []
    for gene in df.loc[df['phase']==phase,'gene_name']:
        #corr_G1.append(counts.loc[gene,].corr(counts_G1))
        corr_S.loc[gene,'S'] = counts.loc[gene,].corr(counts_S)
    corr_S['S'].plot(kind='density')
    plt.show()
    corr_S_filter = corr_S.loc[corr_S['S']>=0.3,]
    genelist_S_filter = pd.DataFrame()
    genelist_S_filter['gene_name'] = corr_S_filter.index
    genelist_S_filter['phase'] = phase
    
    return genelist_S_filter
    

genelist_G1_S_filter = genelist_filter('G1/S')
genelist_G1_S_filter

genelist_S_filter = genelist_filter('S')
genelist_S_filter

genelist_G2_M_filter = genelist_filter('G2/M')
genelist_G2_M_filter

genelist_M_filter = genelist_filter('M')
genelist_M_filter

genelist_M_G1_filter = genelist_filter('M/G1')
genelist_M_G1_filter

genelist_filter = pd.concat([genelist_G1_S_filter,genelist_S_filter,genelist_G2_M_filter,genelist_M_filter,genelist_M_G1_filter],axis=0)
genelist_filter

df_G1_S = pd.merge(genelist_filter[genelist_filter['phase']=='G1/S'],counts,left_on='gene_name',right_on=counts.index,how='left')
df_G1_S = df_G1_S.iloc[:,2:]
df_G1_S.loc['mean',] = df_G1_S.mean(axis=0)
df_G1_S = df_G1_S.loc['mean',]
data_G1_S = pd.DataFrame(df_G1_S)
data_G1_S.columns = ['G1/S']
data_G1_S

df_G1_S = pd.merge(genelist_filter[genelist_filter['phase']=='S'],counts,left_on='gene_name',right_on=counts.index,how='left')
df_G1_S = df_G1_S.iloc[:,2:]
df_G1_S.loc['mean',] = df_G1_S.mean(axis=0)
df_G1_S = df_G1_S.loc['mean',]
data_S = pd.DataFrame(df_G1_S)
data_S.columns = ['S']
data_S

df_G1_S = pd.merge(genelist_filter[genelist_filter['phase']=='G2/M'],counts,left_on='gene_name',right_on=counts.index,how='left')
df_G1_S = df_G1_S.iloc[:,2:]
df_G1_S.loc['mean',] = df_G1_S.mean(axis=0)
df_G1_S = df_G1_S.loc['mean',]
data_G2_M = pd.DataFrame(df_G1_S)
data_G2_M.columns = ['G2/M']
data_G2_M

df_G1_S = pd.merge(genelist_filter[genelist_filter['phase']=='M'],counts,left_on='gene_name',right_on=counts.index,how='left')
df_G1_S = df_G1_S.iloc[:,2:]
df_G1_S.loc['mean',] = df_G1_S.mean(axis=0)
df_G1_S = df_G1_S.loc['mean',]
data_M = pd.DataFrame(df_G1_S)
data_M.columns = ['M']
data_M

df_G1_S = pd.merge(genelist_filter[genelist_filter['phase']=='M/G1'],counts,left_on='gene_name',right_on=counts.index,how='left')
df_G1_S = df_G1_S.iloc[:,2:]
df_G1_S.loc['mean',] = df_G1_S.mean(axis=0)
df_G1_S = df_G1_S.loc['mean',]
data_M_G1 = pd.DataFrame(df_G1_S)
data_M_G1.columns = ['M/G1']
data_M_G1

data = pd.concat([data_G1_S,data_S,data_G2_M,data_M,data_M_G1],axis=1)
data

from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
data[['G1/S', 'S', 'G2/M','M','M/G1']] = scaler.fit_transform(data[['G1/S', 'S', 'G2/M','M','M/G1']])

data['max'] = data.max(axis=1)
data['phase'] = data.iloc[:,0:5].idxmax(axis=1)
data.loc[data['max']<0,'phase'] = 'G0'
data

adata1.obs['cell_cycle_phase'] = data['phase']
adata1.obs.to_csv('~/cellcycle_phase.csv')

