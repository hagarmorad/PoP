#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 09:46:29 2023

@author: hagar

This script was developed with Hagar Lester and Neta to generate reports for nOPV2.
nOPV2 is different from the other poliovirus thus required different analysis.
The output:
    1. mutation table with N's - if you see alot of Ns in the 5'UTR you should suspect it NOT nOPV2, as its 5UTR is longer than other polioviruses'
    2. mutations table with no Ns - in this excle file you can find the nutations with no Ns and the "attenuation" sheet. 
                                    this sheet summarize all the decreases / increases mutations found in this run. updated to Jan-23
    3. nOPV mutations attenuations - this file is derived from the mutations.xlsx file. this is Lester's format to analyze nOPV2 and we automated it for him. 
    
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
    get all mutations in a specific range    

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
    
    return [x for x in mut if start < int(x.split("-" + mut_type)[0][1:-1]) and end > int(x.split("-" + mut_type)[0][1:-1])]


def remove_fp_muts(mut_list, mut_str):
    '''
    ** Help function for lester_format function. **
    *** you should understand lester_format function before reading this one.***
    
    To avoid false positive mutations we want to remove consequetive deletions from the mutations lists.
    we assume these are false as they probably exist:
    1- because someone cut the sequences manually before aligning it
    2- because of bad sequencing quelity.

    Parameters
    ----------
    mut_list : list 
        a mutation list. each item format is: nuc pos nuc - mut type. 
        exmple: G7466A-syn
    mur_str : string
        each mutation in the list contains the mutation type "-syn" or "-nonsyn"
        

    Returns
    -------
    None.

    '''
    counter = 0
    prev = 0
    new_mut_list = []
    for item in mut_list:
        new_mut_list.append(item) # first add the mutation - delete it later if necessery
        mut = item.split(mut_str)[0]
        pos = int(mut[1:-1])
        
        if mut[-1] == "-":
            if prev + 1 == pos: # if this positions follows the previous one
                counter += 1
                if counter >= 5:
                    #remove the current mutation and all pervious 5 from the list.
                    new_mut_list = [x for x in new_mut_list if not int(x.split(mut_str)[0][1:-1]) in [pos,prev, prev-1,prev-2,prev-3,prev-4]] 
                     
            else: # if we got here we are not in consequetive sequence so we reset the counter
                counter = 0 

            prev = pos  # remember position to compare it with next position
        
    return new_mut_list
    
def lester_format(muttbl, ref_name):
    '''
    this section is for Lester. 
    it generates the yellow and the green part of Lester's table    
    
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
            nonsyn = [x for x in nonsyn if not str(x) == "nan"]  #remove None, cast to list

            syn = remove_fp_muts(syn, "-syn")
            nonsyn = remove_fp_muts(nonsyn, "-nonsyn")
            
            # yellow part and green part are different peresntation of the mutations according to Lester's format
            
            # yellow part
            
            CRE_mut = get_mut_in_range(syn, 120, 182, "syn")
            gaps_in_cre = [x for x in CRE_mut if x.split('-syn')[0][-1] == '-' ]   
            final_table.loc["Presence of CRE5 in 5 NCR", sample] = True if len(gaps_in_cre) < 10 else False

            domain_V_mut = get_mut_in_range(syn, 529, 596, "syn")
            final_table.loc["Presence of modified Domain V in 5 NCR", sample] = True if len(domain_V_mut) <= 3 else False
            
            mut_814 =[x for x in nonsyn if x.split("-nonsyn")[0][1:-1] == '814']
            final_table.loc["Presence of A at nucleotide 814", sample] = True if mut_814 else False
    
            mut_817 =[x for x in nonsyn if x.split("-nonsyn")[0][1:-1] == '817']
            final_table.loc["Presence of T at nucleotide 817", sample] = True if mut_817 else False
            
            mut_1375 =[x for x in nonsyn if x.split("-nonsyn")[0][1:-1] == '1375']
            final_table.loc["Presence of T at nucleotide 1375", sample] = True if mut_1375 else False
            
            final_table.loc["Presence of Rec1 K38R mutation in 3D coding sequence", sample] = False if get_mut_in_range(nonsyn, 6158, 6160, "nonsyn") else True
            final_table.loc["Presence of HiFi D53N mutation in 3D coding sequence", sample] = False if get_mut_in_range(nonsyn, 6203, 6205, "nonsyn") else True
            
            ko_cre_mut = get_mut_in_range(syn, 4508, 4560, "syn")
            final_table.loc["Presence of KO CRE"] = True if len(ko_cre_mut) <= 3 else False


            #green part 
            
            final_table.loc["Mutations in CRE5", sample] = ("; ").join(CRE_mut)
            final_table.loc["Mutations in Domain II in 5' NCR",sample] = ("; ").join(get_mut_in_range(syn, 185, 223, "syn"))
            final_table.loc["Mutation U459C in Domain IV in 5 NCR",sample] = 'U459C' if 'T459C-syn' in syn else ""
            final_table.loc["Mutations in modified Domain V in 5 NCR",sample] = ("; ").join(domain_V_mut)
            final_table.loc["Mutation at VP1-143",sample] = ("; ").join(get_mut_in_range(nonsyn, 2969, 2971, "nonsyn"))
            final_table.loc["Mutation at VP1-171",sample] = ("; ").join(get_mut_in_range(nonsyn,3056, 3056, "nonsyn"))
            final_table.loc["Mutation at VP1-295",sample] = ("; ").join(get_mut_in_range(nonsyn, 3428, 3430, "nonsyn"))
            final_table.loc["Mutations in 2C KO CRE",sample] = ("; ").join(get_mut_in_range(nonsyn, 4508, 4560, "nonsyn") )
                       
            final_table.loc["syn",sample] =  ("; ").join(syn) 
            final_table.loc["nonSyn",sample] =  ("; ").join(nonsyn) 
            
 

            
           
            


    
    
    return final_table
    
def run(fasta):
    regions = MAIN_SCIPT_DIR + "refs/nOPV_regions.csv"
    reference = MAIN_SCIPT_DIR + "refs/nOPV.fasta"
    create_dirs(["alignment", "reports"])
    mut_file = "reports/mutations.xlsx"
    mafft(fasta, reference, "alignment/all_aligned.fasta")
    
    signatures.run("alignment/all_aligned.fasta", regions, "reports/mutations.xlsx",noN = 1)
    signatures.run("alignment/all_aligned.fasta", regions, "reports/mutations_with_N.xlsx",noN = 0)

    atten = pd.read_csv(MAIN_SCIPT_DIR + "refs/nOPV2_Attenuation_Mutations.csv")
    
    muts = pd.read_excel(mut_file)
    
    atten = atten.merge(muts, how='left', on='nt_position_on_genome').dropna(thresh=3)
    
    writer = pd.ExcelWriter(mut_file, engine="openpyxl", mode="a", if_sheet_exists="overlay")
    atten.to_excel(writer, index=False, sheet_name= "attenuation")
    writer.save()
    
    muttbl = pd.read_excel("reports/mutations.xlsx", "nOPV_regions")
    muttbl_w_n = pd.read_excel("reports/mutations_with_N.xlsx")
    with open(reference) as f:
        ref_name = f.readline().replace(">","").split(" ")[0]
    lester_format(muttbl, ref_name).to_csv(mut_file.replace("mutations.xlsx", "nOPV_mutations_attenuations.csv"))
    
