import pandas as pd
# import matplotlib.pyplot as plt
import numpy as np
# import seaborn as sns
# sns.set_context("paper", font_scale=2.0)
# sns.set_style('whitegrid')
import math
# import scipy
# import ipywidgets as widgets
# from numba import jit

list_of_aa_codes="ACDEFGHIKLMNPQRSTVWY"
aa = 'Phe-Cys-Tyr-Gln-Leu-Asn-Gly-Ile-Pro-Arg-Ala-Asp-Glu-His-Lys-Met-Ser-Thr-Trp-Val'.split('-')
code_of_aa={'G':'Gly','P':'Pro','A':'Ala','V':'Val','L':'Leu','I':'Ile','M':'Met','C':'Cys','F':'Phe','Y':'Tyr','W':'Trp','H':'His','K':'Lys','R':'Arg','Q':'Gln','N':'Asn','E':'Glu','D':'Asp','S':'Ser','T':'Thr'}

def normalise(df):
    return (df-df.min())/(df.max()-df.min())

def integ_lists (x,y):
    # Trapezium Rule
    xx=x
    yy=y
    temp=0
    for i in range(1,len(xx)):
        h=abs(xx[i]-xx[i-1])
        temp+=h*(yy[i]+yy[i-1])/2
    return temp

def integ_df (x,y):
    # Trapezium Rule
    xx=x.dropna()
    yy=y.dropna()
    temp=0
    for i in range(1,len(xx)):
        h=abs(xx[i]-xx[i-1])
        temp+=h*(yy[i]+yy[i-1])/2
    return temp

def amino_acid_spectrum_simulator(AA, T, nu, f_nu,print_T=False):
    #rnd = np.random.default_rng(rand)
    tot_int = pd.DataFrame(f_nu.sum())
    tot_int.columns = ['Sigma']
    lb, ub = 497.683, 1627.830
    weights=f_nu/tot_int.Sigma.max()
    rel_int=tot_int/tot_int.Sigma.max()
    lam=int(T*rel_int.loc[AA])
    N=np.random.poisson(lam)
    while (N==0):
        N=np.random.poisson(lam)
    draw_nu = nu[AA].sample(N, random_state=N,weights=weights[AA], replace=True)
    draw_f_nu=f_nu[AA][draw_nu.index].values
    final_data=pd.DataFrame()
    final_data['Wavenumber'] = draw_nu.values
    final_data['Intensity'] = draw_f_nu
    if print_T==True:
        print('Total photons at T={} for {} is {}.'.format(T, AA, N))
    return final_data.sort_values(by="Wavenumber").reset_index(drop=True)

    
def amino_acid_seq_sim(seq,T,nu,f_nu):
    seq_split=seq.split('-')
    count=0
    prev_max=0
    tot_photons_em=0
    x_final_em=[]
    y_final_em=[]
    for AA in seq_split:
        sim_dat_AA=amino_acid_spectrum_simulator(AA,T,nu,f_nu)
        x_temp=[x+prev_max for x in list(sim_dat_AA['Wavenumber'].values)]
        x_final_em.extend(x_temp)
        y_final_em.extend(list(sim_dat_AA['Intensity'].values))
        if(count!=0):
            prev_max=max(prev_max, max(x_temp))
        count+=1
    final_data=pd.DataFrame()
    final_data['Wavenumber']=x_final_em
    final_data['Intensity']=y_final_em
    return final_data.sort_values(by="Wavenumber").reset_index(drop=True)

def amino_acid_code_seq_sim(seq_,T,nu,f_nu):
    seq=''
    count=0
    for x in seq_:
        if (count!=0):
            seq+='-'
        seq+=code_of_aa[x]
        count+=1
    return amino_acid_seq_sim(seq,T,nu,f_nu)


def amino_acid_seq_sim_actual(seq,T,nu,f_nu):
    seq_split=seq.split('-')
    count=0
    #prev_max=0
    tot_photons_em=0
    x_final_em=[]
    y_final_em=[]
    for AA in seq_split:
        sim_dat_AA=amino_acid_spectrum_simulator(AA,T,nu,f_nu)
        #print(sim_dat_AA)
        #print(AA, prev_max)
        x_temp=[x for x in list(sim_dat_AA['Wavenumber'].values)]
        x_final_em.extend(x_temp)
        y_final_em.extend(list(sim_dat_AA['Intensity'].values))
        #prev_max=max(prev_max, max(x_temp))
        count+=1
    final_data=pd.DataFrame()
    final_data['Wavenumber']=x_final_em
    final_data['Intensity']=y_final_em
    return final_data

def noise_simulator(T,rel,nu,f_nu):
    lam=int(T*rel)
    N=np.random.poisson(lam)
    while (N==0):
        N=np.random.poisson(lam)
    draw_nu = (nu['Ser'].sample(N, random_state=N,weights=None, replace=True)).values
    draw_f_nu=[]
    for x in draw_nu:
        #currently the noise is random (sadly :[)
        random_value=np.random.uniform(20,100,1)
        draw_f_nu.append(random_value[0])
    final_data=pd.DataFrame()
    final_data['Wavenumber'] = draw_nu
    final_data['Intensity'] = draw_f_nu
    return final_data.sort_values(by="Wavenumber").reset_index(drop=True)

def amino_acid_seq_sim_with_noise(seq,T,nu,f_nu):
    seq_split=seq.split('-')
    count=0
    prev_max=0
    tot_photons_em=0
    x_final_em=[]
    y_final_em=[]
    for AA in seq_split:
        if(AA!='Noise'):
            sim_dat_AA=amino_acid_spectrum_simulator(AA,T,nu,f_nu)
        else:
            rel=np.random.uniform(0.5,1)
            sim_dat_AA=noise_simulator(T,rel,nu,f_nu)
        #print(sim_dat_AA)
        #print(AA, prev_max)
        x_temp=[x+prev_max for x in list(sim_dat_AA['Wavenumber'].values)]
        x_final_em.extend(x_temp)
        y_final_em.extend(list(sim_dat_AA['Intensity'].values))
        prev_max=max(prev_max, max(x_temp))
        count+=1
    final_data=pd.DataFrame()
    final_data['Wavenumber']=x_final_em
    final_data['Intensity']=y_final_em
    return final_data.sort_values(by="Wavenumber").reset_index(drop=True)

def amino_acid_seq_sim_with_noise_actual(seq,T,nu,f_nu):
    seq_split=seq.split('-')
    count=0
    #prev_max=0
    tot_photons_em=0
    x_final_em=[]
    y_final_em=[]
    for AA in seq_split:
        if(AA!='Noise'):
            sim_dat_AA=amino_acid_spectrum_simulator(AA,T,nu,f_nu)
        else:
            rel=np.random.uniform(0.5,1)
            sim_dat_AA=noise_simulator(T,rel,nu,f_nu)
        #print(sim_dat_AA)
        #print(AA, prev_max)
        x_temp=[x for x in list(sim_dat_AA['Wavenumber'].values)]
        x_final_em.extend(x_temp)
        y_final_em.extend(list(sim_dat_AA['Intensity'].values))
        #prev_max=max(prev_max, max(x_temp))
        count+=1
    final_data=pd.DataFrame()
    final_data['Wavenumber']=x_final_em
    final_data['Intensity']=y_final_em
    return final_data
    
def amino_acid_code_seq_sim_with_noise_simple(seq__,T,nu,f_nu):
    seq_=list(seq__)
    seq=''
    prop=0.1
    to_repl=np.random.choice(np.arange(len(seq_))), math.ceil(0.1*len(seq_))
    for x in to_repl:
        seq_[x]='!'
    count=0
    for x in seq_:
        if(count!=0):
            seq+="-"
        if(x=="!"):
            seq+="Noise"
        else:
            seq+=code_of_aa[x]
        count+=1
    #print(seq)
    return amino_acid_seq_sim_with_noise_actual(seq,T,nu,f_nu)

def amino_acid_code_seq_sim_wo_noise_simple(seq__,T,nu,f_nu):
    seq_=list(seq__)
    seq=''
#     prop=0.1
#     to_repl=np.random.choice(np.arange(len(seq_))), math.ceil(0.1*len(seq_))
#     for x in to_repl:
#         seq_[x]='!'
    count=0
    for x in seq_:
        if(count!=0):
            seq+="-"
        if(x=="!"):
            seq+="Noise"
        else:
            seq+=code_of_aa[x]
        count+=1
    #print(seq)
    return amino_acid_seq_sim_with_noise_actual(seq,T,nu,f_nu)

def amino_acid_code_seq_sim_with_noise(seq_,T,nu,f_nu,alpha=0.1,beta=1.5,gamma=0.01):
    #treat this like a markov chain that resets after some state
    #so my states are something like n successes and m failures
    #I start at 0 success
    #If i get a success, i go to 1 success and then see what happens from there
    #Ofcourse, the probability of success should decrease then and probability of failure (thus) increases
    #similarly if i fail now i go to 1 failure
    #else i go to 2 success
    #clear what is going on. The only issue is what about how to change the probabilities?
    #Ofcourse the change in probability of success shouldn't be as drastic as change in probability of failure
    #And the probability of success/failure should intrinisically depend on T
    #Have no real data, so can't predict this? :(
    #One way to deal with this is ofcourse take some prior assumptions and hope it is good enough.
    #so say I have something like probability of failure decreases with T, with probability of failure at T=0 being 1, which is reasonable??
    #how should probability decrease?? something like f=e^(-\alpha T) is a good estimate (vaguely exponential decrease would make sense) where \alpha>0
    #and how should consecutive failures/success increase/decrease
    #say I have one success, then the change of failure is f. If I get another success, let's increase this to 
    last_thing=-1
    #1-failure, 0 success, -1 - just started!
    consec_counter=0
    curr_fail_prob=math.exp(alpha*-T)
    seq=''
    count=0
    for x in seq_:
        if (count!=0):
            seq+='-'
        draw_num=np.random.uniform(0,1)
        if(draw_num<=curr_fail_prob):
            if(last_thing!=0):
                consec_counter=0
            last_thing=0
            curr_fail_prob=math.exp(beta*-consec_counter) * math.exp(alpha*-T)
            consec_counter+=1
            seq+='Noise'
        else:
            if(last_thing!=1):
                consec_counter=0
            last_thing=1
            consec_counter+=1
            curr_fail_prob=min(1,curr_fail_prob*math.exp(gamma*consec_counter))
            seq+=code_of_aa[x]
        count+=1
    return amino_acid_seq_sim_with_noise(seq,T,nu,f_nu)
    
def final_data_maker(df_prot_filtered,nu,f_nu):
    count=0
    final_data=pd.DataFrame()
    final_data['Simulated_Data_Pairs_wo_noise']=None
    final_data['Simulated_Data_Pairs_with_noise']=None
    final_data['Sequence']=None
    max_count= df_prot_filtered.Entry.count()
    
    T=50
    #amino_acid_code_seq_sim_with_noise_simple(seq,T)
    for i,row in df_prot_filtered.iterrows():
        if(count>=0):
            seq_curr_=row['Sequence']
            seq_curr=list(seq_curr_)
            omit=[x for x in seq_curr_ if x not in list_of_aa_codes ]
            if(omit):
                continue
            print(count, row['Length'])
            output_wo_noise=amino_acid_code_seq_sim_wo_noise_simple(seq_curr,T,nu,f_nu)
            output_with_noise=amino_acid_code_seq_sim_with_noise_simple(seq_curr,T,nu,f_nu)
            wave_nums_wo_noise=output_wo_noise.Wavenumber.values
            intensities_wo_noise=output_wo_noise.Intensity.values
            wave_nums_with_noise=output_with_noise.Wavenumber.values
            intensities_with_noise=output_with_noise.Intensity.values
            data_with_noise=[[[wave_nums_with_noise[i], intensities_with_noise[i]] for i in range(len(wave_nums_with_noise))]]
            data_wo_noise=[[[wave_nums_wo_noise[i], intensities_wo_noise[i]] for i in range(len(intensities_wo_noise))]]
            data_df=pd.DataFrame()
            data_df["Simulated_Data_Pairs_wo_noise"]=data_wo_noise
            data_df["Simulated_Data_Pairs_with_noise"]=data_with_noise
            data_df["Sequence"]=seq_curr_
            final_data = pd.concat([final_data,data_df],ignore_index=True)
            count+=1
    return final_data