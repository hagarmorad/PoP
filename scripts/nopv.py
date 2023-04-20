#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 09:46:29 2023

@author: hagar
"""

import os
MAIN_SCIPT_DIR = os.path.dirname(__file__)+'/../'
import numpy as np
import pandas as pd
import sys
sys.path.insert(1, MAIN_SCIPT_DIR)

from scripts.mutations import signatures
from utils.utils import create_dirs, mafft


def get_mut_in_range(mut, start, end, mut_type):
    '''
    

    Parameters
    ----------
    mut : list
        a list od synonimus mutations in the format : 'nuc'-'pos'-'nuc'-'-syn' (A123T-syn).
    start : int
        the 1st nuc pos to look for mutations.
    end : int
        the final nuc pos to look for mutations.
    mut_type: str
        syn / nonsyn

    Returns
    -------
    string
        string of all mutations in this range seperated by ;

    '''
    a = "; ".join([x for x in mut if start < int(x.split("-" + mut_type)[0][1:-1]) and end > int(x.split("-" + mut_type)[0][1:-1])])
    return  a
    
def lester_format(muttbl, ref_name):
    '''
    

    Parameters
    ----------
    muttbl : pandas DataFrame
        the mutation table.

    Returns
    -------
    final table - a table in lesters' format with nOPV attenuations.

    '''
    final_table = pd.DataFrame()

    for col_name in muttbl:
        if col_name.endswith("_NT") and not col_name == ref_name + "_NT":
            muttbl[col_name + "_temp"] = np.where(muttbl[ref_name+ "_NT"] ==  muttbl[col_name]  , None ,muttbl[ref_name+ "_NT"] + muttbl["nt_position_on_genome"].astype(str) + muttbl[col_name] )
        
        if col_name.endswith("_AA") and not col_name == ref_name + "_AA":
            muttbl[col_name + "_temp"] = np.where(muttbl[ref_name + "_AA"] ==  muttbl[col_name]  , None ,\
                                                 muttbl[ref_name + "_AA"] + muttbl["aa_position_on_gene"].astype(str) + muttbl[col_name] )
            sample = col_name.split("_AA")[0]
            
            #synonymus mutation = nucleutide got mutation and amino acid not
            syn = np.where( ( muttbl[sample + "_NT_temp"].notnull()) & (muttbl[sample + "_AA_temp"].isnull()) \
                                                     , muttbl[sample + "_NT_temp"] + "-syn" , None)
            syn = syn[syn != np.array(None)].tolist() #remove None, cast to list
                
              
            #non-synonymus mutation = amino acid got mutation 
                
            nonsyn = np.where(muttbl[sample + "_AA_temp"].notnull() , \
                              muttbl[sample + "_NT_temp"] + "-nonsyn " + muttbl["gene_name"] + "(" + muttbl[sample + "_AA_temp"] + ")", None)
            nonsyn = nonsyn[nonsyn != np.array(None)].tolist() #remove None, cast to list
            
            nonsyn = [x for x in nonsyn if not str(x) == "nan"]
    
    
            final_table.loc["Mutations in CRE5", sample] = get_mut_in_range(syn, 120, 182, "syn")
            final_table.loc["Mutations in Domain II in 5' NCR",sample] = get_mut_in_range(syn, 185, 223, "syn")
            final_table.loc["Mutation U459C in Domain IV in 5 NCR",sample] = 'U459C' if 'T459C-syn' in syn else ""
            final_table.loc["Mutations in modified Domain V in 5 NCR",sample] = get_mut_in_range(syn, 529, 596, "syn")
            final_table.loc["Mutation at VP1-143",sample] = get_mut_in_range(nonsyn, 2969, 2971, "nonsyn")
            final_table.loc["Mutation at VP1-171",sample] = get_mut_in_range(nonsyn,3056, 3056, "nonsyn")
            final_table.loc["Mutation at VP1-295",sample] = get_mut_in_range(nonsyn, 3428, 3430, "nonsyn")
            final_table.loc["Mutations in 2C KO CRE",sample] = get_mut_in_range(nonsyn, 4508, 4560, "nonsyn") 
                       
            final_table.loc["syn",sample] =  ("; ").join(syn) 
            final_table.loc["nonSyn",sample] =  ("; ").join(nonsyn) 
    return final_table
    
def run(fasta):
    regions = MAIN_SCIPT_DIR + "refs/nOPV_regions.csv"
    reference = MAIN_SCIPT_DIR + "refs/nOPV.fasta"
    create_dirs(["alignment", "reports"])
    mut_file = "reports/mutations.xlsx"
    mafft(fasta, reference, "alignment/all_aligned.fasta")
    
    signatures.run("alignment/all_aligned.fasta", regions, "reports/mutations.xlsx")

    atten = pd.read_csv(MAIN_SCIPT_DIR + "refs/nOPV2_Attenuation_Mutations.csv")
    
    muts = pd.read_excel(mut_file)
    
    atten = atten.merge(muts, how='left', on='nt_position_on_genome').dropna(thresh=3)
    
    writer = pd.ExcelWriter(mut_file, engine="openpyxl", mode="a", if_sheet_exists="overlay")
    atten.to_excel(writer, index=False, sheet_name= "attenuation")
    writer.save()
    
    muttbl = pd.read_excel("reports/mutations.xlsx", "nOPV_regions")
    with open(reference) as f:
        ref_name = f.readline().replace(">","").split(" ")[0]
    lester_format(muttbl, ref_name).to_csv(mut_file.replace("mutations.xlsx", "nOPV_mutations_attenuations.csv"))
    
