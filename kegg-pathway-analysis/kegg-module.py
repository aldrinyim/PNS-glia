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
from slugify import slugify, Slugify, UniqueSlugify

current_wd = os.getcwd()
print('working directory: ',current_wd)

kegg_map_add = "/Users/aldrinyim/Documents/database/kegg_05092016/mmu_ensembl.list"
kegg_pathway_add = "/Users/aldrinyim/Documents/database/kegg_05092016/mmu_pathway.list"
map_name_add = "/Users/aldrinyim/Documents/database/kegg_05092016/map_title.tab"

parser = argparse.ArgumentParser()
#parser.add_argument("-s", "--scRNA", help=" ESSENTIAL Single cell RNA-seq matrix")
#parser.add_argument("-c", "--cluster", help=" ESSENTIAL Cluster identity from Seurat")
parser.add_argument("-g", "--gtf", help=" ESSENTIAL GTF file")
args = parser.parse_args()

if args.gtf == None:
    parser.print_help()
    sys.exit()

def convert_ensembl_to_genename(tlist):
    return [mm10_gtf.get_gene_name(x) for x in tlist]

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

selected_pathway_ID = tm.get_all_modules()
print(selected_pathway_ID)
#selected_pathway_ID = ['00600']
#print(convert_ensembl_to_genename(tm.get_ensembl_by_pathwayID(selected_pathway_ID[0])))
outp = open('index.txt','w')
outj = open('index_name.txt','w')
for j in selected_pathway_ID:
    print('working on : ',j)
    filename = j+'_pathwaySymbol.csv'
    outp.write(filename+'\n')
    outf = open(filename,"w")
    outf.write("Pathway,Gene,EnsemblID\n")
    ensembl_list = tm.get_ensembl_by_pathwayID(j)
    #print(ensembl_list)
    wvirgin = 1
    for i in ensembl_list:
        #print(' ---- Gene : ',i)
        if mm10_gtf.get_gene_name(i) is not None:
            custom_slugify = Slugify(to_lower=True)
            custom_slugify.separator = '.'
            slug_map_name = custom_slugify(tm.get_map_name(j))
            print('modified map name : ',slug_map_name)
            if wvirgin == 1:
                outj.write(slug_map_name+'1\n')
                wvirgin = 0
            outf.write(slug_map_name+','+mm10_gtf.get_gene_name(i)+','+i+'\n')
    outf.close()
outp.close()
outj.close()
