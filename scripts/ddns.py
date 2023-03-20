#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 12:13:32 2023

@author: hagar
"""
import pandas as pd 
import csv
import subprocess
from utils.utils import create_dirs, mafft, change_header
import os
import pysam
from io import StringIO
from statistics import mean
from scripts.mutations.signatures import run_ddns

def get_barcode_name(barcode_csv):
    
    with open(barcode_csv, mode='r') as infile:
        reader = csv.reader(infile)
        barcodes = dict((rows[0],rows[1]) for rows in reader)
    
    return barcodes
    


def minimap(fastq_path, reference, barcodes, output):
    minimap = "minimap2 -ax map-ont %(ref)s %(fastq_dir)s/*.fastq.gz | \
                samtools view -bS -F 4 | samtools sort > %(output)s.bam"
    minimap = "minimap2 -ax map-ont %(ref)s %(fastq_dir)s/*.fastq.gz | \
                samtools sort > %(output)s.bam"
    split = "bamtools split -in %(bam)s.bam -reference"
    
    for fastq_dir in os.listdir(fastq_path):
        if os.path.isdir(fastq_path+fastq_dir) and fastq_dir in barcodes.keys():
            sample = barcodes[fastq_dir]
            subprocess.call(minimap % dict(ref=reference, fastq_dir=fastq_path+fastq_dir,\
                                           output=output+"BAM/"+sample), shell=True)
            subprocess.call(split % dict(bam=output+"BAM/"+sample), shell=True)
            
            
def depth(output):
    samtools_index = "samtools index %(bam)s"
    depth = "samtools depth -a %(bam)s > %(output)s.txt"
    for bam in os.listdir(output+"BAM/"):
        if "bai" not in bam:
            subprocess.call(samtools_index % dict(bam=output+"BAM/"+bam), shell=True)
            if "REF_Sabin" in bam:
                subprocess.call(depth % dict(bam=output+"BAM/"+bam, \
                                             output=output+"depth/"+bam.split(".bam")[0]), shell=True) 

def cns(output):
    for bam in os.listdir(output+"BAM/"):
        if "REF_Sabin" in bam:
            sample = bam.split(".bam")[0]
            ivar_cns = "samtools mpileup -A %(bam)s.bam | ivar consensus -t 0.6 -m 1 -p %(cns)s.fa"
            ivar_cns5 = "samtools mpileup -A %(bam)s.bam | ivar consensus -t 0.6 -m 5 -p %(cns)s.fa"
            subprocess.call(ivar_cns % dict(bam=output+"BAM/"+sample, cns=output+"CNS/"+sample), shell=True) 
            subprocess.call(ivar_cns5 % dict(bam=output+"BAM/"+sample, cns=output+"CNS5/"+sample), shell=True) 
            os.remove(output+"CNS/"+sample+".qual.txt")
            os.remove(output+"CNS5/"+sample+".qual.txt")
    change_header(output+"CNS5/")
    change_header(output+"CNS/")
    

def align(s1_ref, s2_ref, s3_ref, output):
    subprocess.call("cat " + output+"CNS5/*REF_Sabin1* > " + output + "alignment/sabin1_not_aligned.fasta", shell=True)
    subprocess.call("cat " + output+"CNS5/*REF_Sabin2* > " + output + "alignment/sabin2_not_aligned.fasta", shell=True)
    subprocess.call("cat " + output+"CNS5/*REF_Sabin3* > " + output + "alignment/sabin3_not_aligned.fasta", shell=True)
    mafft(output+"alignment/sabin1_not_aligned.fasta", s1_ref, output+"alignment/sabin1_aligned.fasta")     
    mafft(output+"alignment/sabin2_not_aligned.fasta", s2_ref, output+"alignment/sabin2_aligned.fasta")
    mafft(output+"alignment/sabin3_not_aligned.fasta", s3_ref, output+"alignment/sabin3_aligned.fasta")


def get_stat(depth_file, sample, stat, ref):
    
    mapped_reads = int(stat.loc[ref]["numreads"])
    covered_bases = int(stat.loc[ref]["covbases"])
    coverage = int(stat.loc[ref]["coverage"])
    
    depths = [int(x.split('\t')[2]) for x in open(depth_file).readlines()]
    mean_depth= str(round(mean(depths),2)) if depths else ''
    min_depth = min(depths) if depths else ''
    max_depth = max(depths) if depths else ''
    genome_size = len(depths)
    
    breadth_cns5 = len([i for i in depths if i > 5])
    coverage_cns5 = round(breadth_cns5 / genome_size * 100,2)  if genome_size else ''
    breadth_cns1 = len([i for i in depths if i > 1])
    coverage_cns1 = round(breadth_cns1 / genome_size * 100,2)  if genome_size else ''
    
    return mapped_reads, covered_bases, coverage, mean_depth, min_depth, max_depth, coverage_cns1, coverage_cns5


def QC_reports(barcodes, output):
    general_qc = pd.DataFrame(columns=["total_reads", "mapped_reads", "%mapped", "%Sabin1 (from mapped_reads)", "%Sabin2 (from mapped_reads)", "%Sabin3 (from mapped_reads)"])
    sabin1_df = pd.DataFrame(columns=["mapped_reads", "covered_bases", "coverage", "mean_depth", "min_depth", "max_depth", "coverage", "coverage_cns5"])
    sabin2_df = pd.DataFrame(columns=["mapped_reads", "covered_bases", "coverage", "mean_depth", "min_depth", "max_depth", "coverage", "coverage_cns5"])
    sabin3_df = pd.DataFrame(columns=["mapped_reads", "covered_bases", "coverage", "mean_depth", "min_depth", "max_depth", "coverage", "coverage_cns5"])
    
    for sample in barcodes.values():
        
        #general information
        total_reads = pysam.AlignmentFile(output + "BAM/" + sample + ".bam").count(until_eof=True)
        mapped_reads = total_reads - pysam.AlignmentFile(output + "BAM/" + sample + ".REF_unmapped.bam").count(until_eof=True)
        mapped_percentage = round(mapped_reads / total_reads * 100, 2)
        stat = pd.read_csv(StringIO(pysam.coverage(output + "BAM/" + sample + ".bam")), sep='\t').set_index("#rname")
        sabin1_reads_pre = (stat.loc["Sabin1"]["numreads"]/mapped_reads*100)
        sabin2_reads_pre = (stat.loc["Sabin2"]["numreads"]/mapped_reads*100)
        sabin3_reads_pre = (stat.loc["Sabin3"]["numreads"]/mapped_reads*100)
       
        general_qc.loc[sample] = (total_reads, mapped_reads, mapped_percentage, sabin1_reads_pre, sabin2_reads_pre, sabin3_reads_pre)
        
        # info by reference
        depth_file = output + "depth/" + sample +".REF_Sabin1.txt"
        sabin1_df.loc[sample] = get_stat(depth_file, sample, stat, "Sabin1" )
        depth_file = output + "depth/" + sample +".REF_Sabin2.txt"
        sabin2_df.loc[sample] = get_stat(depth_file, sample, stat, "Sabin2" )
        depth_file = output + "depth/" + sample +".REF_Sabin3.txt"
        sabin3_df.loc[sample] = get_stat(depth_file, sample, stat, "Sabin3" )
        
    with pd.ExcelWriter(output + "reports/" + 'reports.xlsx') as writer:
        general_qc.to_excel(writer, sheet_name='General_QC')
        sabin1_df.to_excel(writer, sheet_name='Sabin1')
        sabin2_df.to_excel(writer, sheet_name='Sabin2')
        sabin3_df.to_excel(writer, sheet_name='Sabin3')
    

def run(fastq_path, barcode_csv, output):
    par_dir = os.path.dirname(os.path.abspath(os.path.dirname(__file__ )))
    reference = par_dir + "/refs/VP1.fasta"
    s1_ref = par_dir + "/refs/VP1_sabin1.fasta"
    s2_ref = par_dir + "/refs/VP1_sabin2.fasta"
    s3_ref = par_dir + "/refs/VP1_sabin3.fasta"
    create_dirs([output+"BAM/", output+"depth/", output+"CNS/", output+"CNS5/", output+"alignment/", output + "reports/"])
    barcodes = get_barcode_name(barcode_csv)
    # minimap(fastq_path, reference, barcodes, output)
    # depth(output)
    # cns(output)
    # align(s1_ref, s2_ref, s3_ref, output)
    # QC_reports(barcodes, output)
    
    run_ddns(output + "alignment/sabin1_aligned.fasta" ,par_dir + "/refs/sabin1.csv", output + "reports/sabin1_mutations.xlsx")
    run_ddns(output + "alignment/sabin2_aligned.fasta" ,par_dir + "/refs/sabin2.csv", output + "reports/sabin2_mutations.xlsx")
    run_ddns(output + "alignment/sabin3_aligned.fasta" ,par_dir + "/refs/sabin3.csv", output + "reports/sabin3_mutations.xlsx")
    

