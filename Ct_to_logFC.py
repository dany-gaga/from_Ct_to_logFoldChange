#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 11:16:07 2023

@author: dukpe
"""

# import libraries
import sys
import pandas as pd
import numpy as np
import math as m

# read given an csv file argument as file
file = pd.read_csv(sys.argv[1], sep=',') # sep could be change in case

# create three new columns
# Each new column contain a corrected value for the genes Cts, the housekepping genes 1 and 2 Cts 
# using their respctive efficiency (efficieency corrected Ct)

file['cor_Gene'], file['cor_HK1'], file['cor_HK2'] = [(file.iloc[:,3]*((np.log(file.iloc[:,6]+1)/np.log(2)))),
                                                  (file.iloc[:,4]*((np.log(file.iloc[:,7]+1)/np.log(2)))),
                                                  (file.iloc[:,5]*((np.log(file.iloc[:,8]+1)/np.log(2))))]

# Here the mean for each technical replicates group is calculated
cor_ct_mean = file.groupby([file.iloc[:,0],file.iloc[:,1], file.iloc[:,2]]).aggregate(
    {'cor_Gene':'mean', 'cor_HK1': 'mean', 'cor_HK2':'mean'})

# This part generate the deltaCt using this formula: dCt (target_Ct - Housekeeping_Ct)
cor_ct_mean['gene_HK1'], cor_ct_mean['gene_HK2'] = [(cor_ct_mean['cor_Gene']-cor_ct_mean['cor_HK1']),
                                           (cor_ct_mean['cor_Gene']-cor_ct_mean['cor_HK2'])]

# deltaCt (mean value for >1 control gene)
cor_ct_mean['dCt']=((cor_ct_mean['gene_HK1'] + cor_ct_mean['gene_HK2']) / 2)
cor_ct_mean = cor_ct_mean.reset_index()

# Mean for the three biological replicate
dCt_mean = cor_ct_mean.groupby([cor_ct_mean.iloc[:,0],cor_ct_mean.iloc[:,1]]).aggregate({'dCt':'mean'})

# save the dCt as a csv file
dCt_mean.to_csv('dcT_mean.csv')
dCt_m = dCt_mean.reset_index()

gene_list = list(set(list(dCt_m["Genes"].to_list()))) # takes all genes present in the study and store it as a list

logFC_table = pd.DataFrame(gene_list, columns=['Genes']) # creating a onecolumn dataframe with the list obtained above

# the logFC_table_construction function help into caluclating the logFC for different comparison using 
# this formula : log2(2^-((resistant dCt) - (sensisble or control dCT)))
def logFC_table_construction(ag, tbl, r, c, **kwargs):
    s = kwargs.get('s', None)
 
    if s != None:
        r_s = []
        for value in tbl.iloc[:,0]:
            if value == value:
                r_s.append(np.log((2**-((ag.loc[value,r]['dCt']) - (ag.loc[value,s]['dCt']))),2))

        c_s = []
        for value in tbl.iloc[:,0]:
            if value == value:
                c_s.append(np.log((2**-((ag.loc[value,c]['dCt']) - (ag.loc[value,s]['dCt']))),2))
         
        r_c = []
        for value in tbl.iloc[:,0]:
            if value == value:
                r_c.append(np.log((2**-((ag.loc[value,r]['dCt']) - (ag.loc[value,c]['dCt']))),2))
        
        tbl["logFC_r_c"], tbl["logFC_r_s"], tbl["logFC_c_s"] = r_c, r_s, c_s
        return(tbl)
    
    elif type(r) == list:
        for i in r:
            vn = "c_" + str(i)
            vn = []
            for value in tbl.iloc[:,0]:
                if value == value:
                    vn.append(np.log((2**-((ag.loc[value,i]['dCt']) - (ag.loc[value,c]['dCt']))),2))
            tbl['logFC_'+str(i)] = vn
        return(tbl)
    else:
        r_c = []
        for value in tbl.iloc[:,0]:
            if value == value:
                r_c.append(np.log((2**-((ag.loc[value,r]['dCt']) - (ag.loc[value,c]['dCt']))),2))
        tbl["logFC_r_c"] = r_c
        return(tbl)

# here the user is asked for the sets he need the FC within his data
# here the user is asked for the sets he need the FC within his data
quest = input('do you have multiple  resistant group yes or no: ')

if quest == 'no':
    
    R = str(input("What is your resitant group name: "))
    C = str(input("What is your control group name: "))
    S = str(input("What is your susceptible group name(write na if unexistant): "))
    
else:
    R = list(map(str, input("What are your resitant group names: ").split()))
    C = str(input("What is your control group name: "))
    S = str(input("What is your susceptible group name(write na if unexistant): "))

if S == 'na':
    t = logFC_table_construction(dCt_mean, logFC_table, r= R, c= C) 

else:
    t = logFC_table_construction(dCt_mean, logFC_table, r= R, c= C, s=S)
print(t)
# le logFC is generated as a csv file
t.to_csv('lfc_from_Ct.csv')
