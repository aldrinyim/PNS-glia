#!/usr/bin/python3

'''
Private Class, do not modify
'''
from collections import defaultdict
import numpy as np

class Gtf:

    def __init__(self):
        self.gene_general_info = defaultdict(list)
        self.genes_exons = defaultdict(dict)
        self.genes_introns = defaultdict(dict)
        self.genes_refined_exons = defaultdict(dict)
        self.gene_ids = []

    def zero_runs(self,a,value):
        # Create an array that is 1 where a is 0, and pad each end with an extra 0.
        iszero = np.concatenate(([0], np.equal(a, value).view(np.int8), [0]))
        absdiff = np.abs(np.diff(iszero))
        # Runs start and end where absdiff is 1.
        ranges = np.where(absdiff == 1)[0].reshape(-1, 2).tolist()
        for i in range(0,len(ranges)):
            ranges[i] = [ranges[i][0],ranges[i][1]-1]
        return ranges

    def create_intron_map(self):
        for gene in self.gene_ids:
            #pull out all the exons for one gene
            specific_gene_exons_id = self.genes_exons[gene].keys()
            strand = self.gene_general_info[gene][5]
            #position w.r.t. to genome
            specific_gene_start = self.gene_general_info[gene][3]
            specific_gene_end = self.gene_general_info[gene][4]
            #Create list of 0
            listofzeros = [0] * (specific_gene_end-specific_gene_start+1)
            for exon1 in specific_gene_exons_id:                             # Zero-based coordinates, no need to add 1
                es = self.genes_exons[gene][exon1][2] - specific_gene_start
                en = self.genes_exons[gene][exon1][3] - specific_gene_start
                while es <= en:
                    listofzeros[es]+=1
                    es+=1
            #Extracting zero regions
            range_of_introns = self.zero_runs(np.asarray(listofzeros),0)
            #Extracting introns
            for i in range(0,len(range_of_introns)):
                if i == 0:
                    true_ig_s = specific_gene_start + range_of_introns[i][0]
                    true_ig_e = specific_gene_start + range_of_introns[i][1]
                    self.genes_introns[gene]=[[true_ig_s,true_ig_e]]
                else:
                    true_ig_s = specific_gene_start + range_of_introns[i][0]
                    true_ig_e = specific_gene_start + range_of_introns[i][1]
                    self.genes_introns[gene].append([true_ig_s,true_ig_e])

    def create_exon_map(self):
        for gene in self.gene_ids:
            #pull out all the exons for one gene
            specific_gene_exons_id = self.genes_exons[gene].keys()
            strand = self.gene_general_info[gene][5]
            #position w.r.t. to genome
            #specific_gene_name = self.gene_general_info[gene][0]
            #specific_gene_chrom = self.gene_general_info[gene][2]
            #specific_gene_strand = self.gene_general_info[gene][5]
            specific_gene_start = self.gene_general_info[gene][3]
            specific_gene_end = self.gene_general_info[gene][4]
            #Create list of 0
            listofone = [0] * (specific_gene_end-specific_gene_start+1)
            for exon1 in specific_gene_exons_id:                             # Zero-based coordinates, no need to add 1
                es = self.genes_exons[gene][exon1][2] - specific_gene_start
                en = self.genes_exons[gene][exon1][3] - specific_gene_start
                while es <= en:
                    if listofone[es] == 0:
                        listofone[es]+=1
                    es+=1
            #Extracting ONE regions
            range_of_exons = self.zero_runs(np.asarray(listofone),1)
            #Extracting introns
            for i in range(0,len(range_of_exons)):
                if i == 0:
                    true_ig_s = specific_gene_start + range_of_exons[i][0]
                    true_ig_e = specific_gene_start + range_of_exons[i][1]
                    self.genes_refined_exons[gene]=[[true_ig_s,true_ig_e]]
                else:
                    true_ig_s = specific_gene_start + range_of_exons[i][0]
                    true_ig_e = specific_gene_start + range_of_exons[i][1]
                    self.genes_refined_exons[gene].append([true_ig_s,true_ig_e])

    def get_intron_map(self):
        return self.genes_introns

    def get_exon_map(self):
        return self.genes_refined_exons

    def get_intron_map_specific_gene(self,gene):
        return self.genes_introns[gene]

    def get_exon_map_specific_gene(self,gene):
        return self.genes_refined_exons[gene]

    def put_gene(self,datain):
        #Extract info from data
        data = datain.split('\t')
        chrom = data[0]
        gtype = data[2]
        start = int(data[3])-1  #Pysam coordinate - zero-based
        end = int(data[4])-1   #Pysam coordinate - zero-based
        strand = data[6]
        detailed_info = data[8].split('; ')
        #print(data[8])
        gene_id = ''
        gene_name = ''
        gene_biotype = ''
        exon_number = -1
        exon_id = ''

        for seg in detailed_info:
            seg_list = seg.split(' ')
            #print(seg_list)
            if seg_list[0] == "gene_id":
                gene_id = self.remove_quote(seg_list[1]).split('.')[0]
            if seg_list[0] == "gene_name":
                gene_name = self.remove_quote(seg_list[1])
            if seg_list[0] == "gene_biotype":
                gene_biotype = self.remove_quote(seg_list[1])
            if seg_list[0] == "exon_number":
                exon_number = seg_list[1]
            if seg_list[0] == "exon_id":
                exon_id = self.remove_quote(seg_list[1])

        if gtype == 'gene':
            self.gene_general_info[gene_id]=[gene_name,gene_biotype,chrom,int(start),int(end),strand]
            self.gene_ids.append(gene_id)

        elif gtype == 'exon':
            self.genes_exons[gene_id][exon_id]=[exon_number,chrom,int(start),int(end),strand]

    def get_total_ids(self):
        return self.gene_ids

    def get_gene_info(self,gene_id):
        return self.gene_general_info[gene_id]

    def get_gene_no_of_exons(self,gene_id):
        return len(self.genes_exons[gene_id])

    def get_gene_chrom(self,gene_id):
        return self.gene_general_info[gene_id][2]
    def get_gene_start(self,gene_id):
        return self.gene_general_info[gene_id][3]
    def get_gene_end(self,gene_id):
        return self.gene_general_info[gene_id][4]
    def get_gene_strand(self,gene_id):
        return self.gene_general_info[gene_id][5]
    def get_gene_name(self,gene_id):
        if gene_id in self.gene_general_info.keys():
            return self.gene_general_info[gene_id][0]
        else:
            return None
    def remove_quote(self,word):
        return word.split('\"')[1]
