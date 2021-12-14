#!/usr/bin/python3
"""
Author: Aldrin Yim
Version 0.1 (alpha)
Date: Oct 2019
Usage: python3 schwann_expression_visualization.py --help

"""
import matplotlib
matplotlib.use('Agg')
import sys, os, subprocess, argparse, pysam, pandas as pd, seaborn as sns, matplotlib.pyplot as plt, math, statistics, numpy as np
from collections import ChainMap, defaultdict
from itertools import combinations
import multiprocessing, tqdm, itertools, csv

current_wd = os.getcwd()
print('working directory: ',current_wd)

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sndeg", help=" ESSENTIAL SN DEG matrix from edgeR")
parser.add_argument("-v", "--vgdeg", help=" ESSENTIAL VG DEG matrix from edgeR")
parser.add_argument("-n", "--normalized", help=" ESSENTIAL TMM matrix or any normalized matrix")
parser.add_argument("-f", "--snfc", help=" ESSENTIAL SN Fold changes from edgeR")
parser.add_argument("-g", "--vgfc", help=" ESSENTIAL VG Fold changes from edgeR")

args = parser.parse_args()

if args.normalized == None:
    parser.print_help()
    sys.exit()

selected_genes = []

all_coors_x = []
all_coors_y = []

sn_coors_x = []
sn_coors_y = []
vg_coors_x = []
vg_coors_y = []
coup_coors_x = []
coup_coors_y = []

init=0
with open(args.sndeg) as f:
    for line in f:
        readin = line.rstrip()
        data = readin.split('\t')
        if init == 0:
            init = 1
        else:
            gene = data[0]
            selected_genes.append(gene)
f.close()

init=0
with open(args.vgdeg) as f:
    for line in f:
        readin = line.rstrip()
        data = readin.split('\t')
        if init == 0:
            init = 1
        else:
            gene = data[0]
            selected_genes.append(gene)
f.close()

init=0
header = []
sn_genes = []
vg_genes = []

with open(args.snfc) as f:
    for line in f:
        readin = line.rstrip()
        data = readin.split('\t')
        if init == 0:
            header = data
            init=1
        else:
            gene = data[0]
            pvalue = float(data[5])
            logCPM = float(data[4])
            fdr = float(data[6])
            logfc = float(data[3])
            if logfc <= -1 and pvalue <= 0.0001 and logCPM >= 3:
            #if logfc >= 2 and pvalue <= 0.0001 and logCPM >= 4:
                sn_genes.append(gene)
f.close()

init=0
with open(args.vgfc) as f:
    for line in f:
        readin = line.rstrip()
        data = readin.split('\t')
        if init == 0:
            header = data
            init=1
        else:
            gene = data[0]
            pvalue = float(data[5])
            logCPM = float(data[4])
            fdr = float(data[6])
            logfc = float(data[3])
            if logfc >= 1 and pvalue <= 0.0001 and logCPM >= 3:
                vg_genes.append(gene)
                #outf.write('VG'+'|'+gene+'|'+str(pvalue)+'|'+str(logCPM)+'|'+str(logfc)+'\n')
f.close()

coup_genes = list(set(sn_genes).intersection(vg_genes))
sn_genes = [x for x in sn_genes if x not in coup_genes]
vg_genes = [x for x in vg_genes if x not in coup_genes]

print('# of co-up genes : ', len(coup_genes))
print('# of SN genes : ', len(sn_genes))
print('# of VG genes : ', len(vg_genes))

#outf = open("selected_genes.txt",'w')
#for i in coup_genes:
#    outf.write('Co-up'+'\t'+i+'\n')
#for i in vg_genes:
#    outf.write('VG'+'\t'+i+'\n')
#for i in sn_genes:
#    outf.write('SN'+'\t'+i+'\n')
#outf.close()

init=0
header = []
#annotated_genes = ['Mpz','Entpd2',"Gfap",'Pmp2','Sparc','S100a6','S100b','Ncamp','Scn7a','Pdgfra',"Cx3cr1","Acta2"]
#annotated_genes = ['S100b','Sox2','S100a6','Cadm2','Sox10',"Ncam1",'Ldhb','Pmp22','Sparc','Adam10','Pdgfra','Cx3cr1','Acta2','Scn7a','Entpd2','Gfap','Ngfr','Acta2']
#annotated_genes = ['Foxd3','Sox10','Kcna2','Cspg4','Adam10','Mbp','S100a10','Sox2','Mpz','Cldn19','Pdgfrb','Tmem119','Mmrn1','Mecom','Des']
annotated_genes = ['Cldn14','Sox10','Kcna2','Cspg4','Adam10','Mbp','S100a10','Sox2','Mpz','Cldn19','Pdgfrb','Tmem119','Mmrn1','Mecom','Des']

annotated_order = []
annotated_genes_coors_x = []
annotated_genes_coors_y = []
outf = open("selected_genes_enrichemnt-20Nov.txt",'w')
outf.write('Type\tGene\tSN_ratio\tVG_ratio\n')
with open(args.normalized) as f:
    for line in f:
        readin = line.rstrip()
        data = readin.split('\t')
        if init == 0:
            header = data
            init=1
        else:
            gene = data[0]
            trim_g = gene.split('|')[2]

            #sn_whole = (float(data[1])+float(data[2])+float(data[3]))/3+1
            sn_desheathed = (float(data[4])+float(data[5])+float(data[6]))/3+1
            sn_ribo = (float(data[7])+float(data[8])+float(data[9]))/3+1
            vg_whole = (float(data[13])+float(data[14])+float(data[15]))/3+1
            vg_ribo = (float(data[10])+float(data[11])+float(data[12]))/3+1

            if (sn_desheathed >= 0 or sn_ribo >= 0) and (vg_whole >= 0 or vg_ribo >= 0):

                #sn_ratio = np.log2(sn_ribo/sn_whole)
                sn_ratio = np.log2(sn_ribo/sn_desheathed)
                vg_ratio = np.log2(vg_ribo/vg_whole)

                if trim_g == 'Cldn14':
                    print(trim_g,"|",sn_ratio,"|",vg_ratio)

                if gene in selected_genes:
                    if trim_g in annotated_genes:
                        annotated_order.append(trim_g)
                        annotated_genes_coors_x.append(sn_ratio)
                        annotated_genes_coors_y.append(vg_ratio)
                    if gene in sn_genes:
                        sn_coors_x.append(sn_ratio)
                        sn_coors_y.append(vg_ratio)
                        outf.write('SN'+'\t'+gene+'\t'+str(sn_ratio)+'\t'+str(vg_ratio)+'\n')
                    elif gene in vg_genes:
                        vg_coors_x.append(sn_ratio)
                        vg_coors_y.append(vg_ratio)
                        outf.write('VG'+'\t'+gene+'\t'+str(sn_ratio)+'\t'+str(vg_ratio)+'\n')
                    elif gene in coup_genes:
                        coup_coors_x.append(sn_ratio)
                        coup_coors_y.append(vg_ratio)
                        outf.write('Co-UP'+'\t'+gene+'\t'+str(sn_ratio)+'\t'+str(vg_ratio)+'\n')
                    else:
                        outf.write('Others'+'\t'+gene+'\t'+str(sn_ratio)+'\t'+str(vg_ratio)+'\n')
                    all_coors_x.append(sn_ratio)
                    all_coors_y.append(vg_ratio)
f.close()
outf.close()
x = all_coors_x
y = all_coors_y

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
ax.axis([math.floor(min(x))-0.5,math.ceil(max(x))+0.5,math.floor(min(y))-0.5,math.ceil(max(y))+0.5])
dot_size = 15

ax.grid(False)
ax.axhline(y=0,alpha=0.8,color='k',linewidth=1)
ax.axvline(x=0,alpha=0.8,color='k',linewidth=1)
ax.set_facecolor('white')
plt.scatter(all_coors_x,all_coors_y, s=dot_size, alpha=0.2, c='black')
plt.scatter(coup_coors_x,coup_coors_y, s=dot_size, alpha=0.5, c='green')
plt.scatter(sn_coors_x,sn_coors_y, s=dot_size, alpha=0.4, c='red')
plt.scatter(vg_coors_x,vg_coors_y, s=dot_size, alpha=0.4, c='blue')
for i in range(0,len(annotated_order)):
    ax.annotate(annotated_order[i],(annotated_genes_coors_x[i],annotated_genes_coors_y[i]),fontsize=12)

plt.savefig('sciatic_expressivity.pdf')
