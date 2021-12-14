#!/usr/bin/python3
"""
Author: Aldrin Yim
Version 0.1
Date: Oct 2020
Usage: python3 module-based_analysis.py --help
"""
import sys, os, subprocess, argparse, pandas as pd, seaborn as sns, matplotlib.pyplot as plt, math, statistics, numpy as np
from collections import ChainMap, defaultdict
import multiprocessing, tqdm
from gtf import Gtf
from transcriptional_modules import Transcriptional_modules
from scProfiles import ScProfiles
from operator import itemgetter
from scipy import stats

current_wd = os.getcwd()
print('working directory: ',current_wd)

kegg_map_add = "/Users/aldrinyim/Documents/database/kegg_05092016/mmu_ensembl.list"
kegg_pathway_add = "/Users/aldrinyim/Documents/database/kegg_05092016/mmu_pathway.list"
map_name_add = "/Users/aldrinyim/Documents/database/kegg_05092016/map_title.tab"

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--scRNA", help=" ESSENTIAL Single cell RNA-seq matrix")
parser.add_argument("-c", "--cluster", help=" ESSENTIAL Cluster identity from Seurat")
parser.add_argument("-g", "--gtf", help=" ESSENTIAL GTF file")
args = parser.parse_args()

if args.scRNA == None or args.gtf == None:
    parser.print_help()
    sys.exit()

#def convert_ensembl_to_genename(tlist):
#    return [mm10_gtf.get_gene_name(x) for x in tlist]

mm10_gtf = Gtf()
with open(args.gtf) as f:
    for line in f:
        readin=line.rstrip()
        if readin[0] != '#':
            #print(readin)
            mm10_gtf.put_gene(readin)
f.close()

tm = Transcriptional_modules()
tm.put_entrez_ensembl(kegg_map_add)
tm.put_entrez_pathway(kegg_pathway_add)
tm.put_map_name(map_name_add)

#print(tm.get_gene_by_pathwayID('00010'))
#print(tm.get_ensembl_by_pathwayID('00010'))
#print(convert_ensembl_to_genename(tm.get_ensembl_by_pathwayID('00010')))

#Read in scRNA-seq expression matrix and classify them into different Transcriptional_modules
scp = ScProfiles()
header = True
with open(args.scRNA) as f:
    for line in f:
        data = line.rstrip().split(',')
        if header is True:
            scp.define_cells(data)
            header = False
        else:
            scp.put_gene_for_cells(data)
            if data[0] == "Prx":
                #Answer = 2.68903190262809
                print("Prx expression : ", data[1])
f.close()

scp.define_cell_modules(tm,mm10_gtf)
print("------Completed-------")
#print("Mpz pathway")
#print(scp.get_cell_expression_profile_by_moduleID('SN1_AAGTCTGTCAGTCCCT','04514'))
#print("Scn7a pathway")
#print(scp.get_cell_expression_profile_by_moduleID('SN1_AAGTCTGTCAGTCCCT','04261'))

#Put in cluster information
header = True
with open(args.cluster) as f:
    for line in f:
        data = line.rstrip().split(',')
        if header == True:
            header = False
        else:
            scp.put_cluster_info(data)
f.close()

print(scp.get_expression_modules_by_cluster(12,'04261'))
