from functions_for_cluster import *
import math
import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# sns.set_context("paper", font_scale=2.0)
# sns.set_style('whitegrid')
# import scipy
# import ipywidgets as widgets
# from numba import jit

list_of_aa_codes="ACDEFGHIKLMNPQRSTVWY"
aa = 'Phe-Cys-Tyr-Gln-Leu-Asn-Gly-Ile-Pro-Arg-Ala-Asp-Glu-His-Lys-Met-Ser-Thr-Trp-Val'.split('-')
code_of_aa={'G':'Gly','P':'Pro','A':'Ala','V':'Val','L':'Leu','I':'Ile','M':'Met','C':'Cys','F':'Phe','Y':'Tyr','W':'Trp','H':'His','K':'Lys','R':'Arg','Q':'Gln','N':'Asn','E':'Glu','D':'Asp','S':'Ser','T':'Thr'}


df = pd.read_csv('data/Raman_Spectra_20_AA.csv')
aa = [i for i in df.columns if 'Unnamed' not in i]
sfts = [i for i in df.columns if i not in aa]

# Wavenumbers
nu_ = df[sfts].copy()
nu_.columns = aa
# Intensity - a function of wave number
f_nu_ = df[aa].copy()

f_nu = f_nu_.copy()
nu = nu_.copy()

tot_int = pd.DataFrame(f_nu.sum())
tot_int.columns = ['Sigma']

lb, ub = 497.683, 1627.830
weights=f_nu/tot_int.Sigma.max()
rel_int=tot_int/tot_int.Sigma.max()

df_prot = pd.read_csv('data/uniprot_1.tsv',sep='\t')
df_prot.drop('Reviewed',axis=1,inplace=True)
df_prot.drop('Entry Name',axis=1,inplace=True)
df_prot.drop('Protein names',axis=1,inplace=True)
df_prot.drop('Gene Names',axis=1,inplace=True)
df_prot.drop('Organism',axis=1,inplace=True)
df_prot=df_prot.sort_values('Length')

q_hi  = df_prot["Length"].quantile(0.8)
df_prot_filtered = df_prot[(df_prot["Length"] < q_hi)]
df_prot_filtered=df_prot_filtered.sort_values(['Length'], ascending=[False])

final_data=final_data_maker(df_prot_filtered,nu,f_nu)

final_data.to_csv('sim_data_uniprot_seq.csv')
final_data.to_pickle('sim_data_uniprot.pkl.gz')

