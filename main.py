#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 08:51:09 2023

@author: hagar
"""

from utils import parse_input
from  scripts import ddns, nopv


args = parse_input.run()

if args.ddns:
    ddns.run(args.input, args.barcodes, args.output_dir)
if args.nopv:
    # nopv.run("/mnt/project1/projects/POLIO/VDPV2_2022/run29/nOPV_pipe/alignment/run29_nOPV.fasta", "../refs/nOPV_regions.csv", "../refs/nOPV.fasta")
    nopv.run(args.input)