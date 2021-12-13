#!/usr/bin/python3

'''
Private Class, do not modify
'''
from collections import defaultdict

class Transcriptome_profile:

    def __init__(self):
        self.direct_deter_full_expression_profile = {} #
        self.direct_genes = set()
        self.direct_cell_barcode = set()
        self.full_expression_profile = defaultdict(dict)
        self.exon_expression_profile = {}
        self.intron_expression_profile = {}
        self.exon_genes = set()
        self.exon_cell_barcode = set()
        self.intron_genes = set()
        self.intron_cell_barcode = set()
        self.all_cell_barcode = set()
        self.all_genes = set()

    def put_full_gene_profile(self,profile):
        self.direct_deter_full_expression_profile = profile
        for gene in self.direct_deter_full_expression_profile.keys():
            self.direct_genes.add(gene)
            for cell in self.direct_deter_full_expression_profile[gene].keys():
                self.direct_cell_barcode.add(cell)

    def put_exon_expression_profile(self,profile):
        self.exon_expression_profile = profile
        for gene in self.exon_expression_profile.keys():
            self.exon_genes.add(gene)
            self.all_genes.add(gene)
            for cell in self.exon_expression_profile[gene].keys():
                self.all_cell_barcode.add(cell)
                self.exon_cell_barcode.add(cell)
                #print('in exon : ', gene,'|',cell,'|',self.exon_expression_profile[gene][cell])

    def put_intron_expresion_profile(self,inprofile):
        self.intron_expression_profile = inprofile
        for gene in self.intron_expression_profile.keys():
            self.all_genes.add(gene)
            self.intron_genes.add(gene)
            for cell in self.intron_expression_profile[gene].keys():
                self.all_cell_barcode.add(cell)
                self.intron_cell_barcode.add(cell)
                #print('in intron : ', gene,'|',cell,'|',self.intron_expression_profile[gene][cell])

    def combine_expression_profile(self):
        for gene in self.all_genes:
            for barcode in self.all_cell_barcode:
                #print(gene,'|',barcode,'|',self.has_specific_gene_cell_expression_with_feature(gene,barcode,'exon'),'|',self.has_specific_gene_cell_expression_with_feature(gene,barcode,'intron'))
                #print(barcode,'|',self.exon_expression_profile[gene][barcode],'|',self.intron_expression_profile[gene][barcode])
                if self.has_specific_gene_cell_expression_with_feature(gene,barcode,'exon') and self.has_specific_gene_cell_expression_with_feature(gene,barcode,'intron'):
                    if barcode not in self.full_expression_profile[gene]:
                        #self.full_expression_profile[gene][barcode] = set()
                        self.full_expression_profile[gene][barcode] = self.exon_expression_profile[gene][barcode].union(self.intron_expression_profile[gene][barcode])
                    else:
                        self.full_expression_profile[gene][barcode] = self.full_expression_profile[gene][barcode].union(self.exon_expression_profile[gene][barcode].union(self.intron_expression_profile[gene][barcode]))
                    #self.full_expression_profile[gene][barcode] = self.exon_expression_profile[gene][barcode]+self.intron_expression_profile[gene][barcode]
                elif self.has_specific_gene_cell_expression_with_feature(gene,barcode,'exon'):
                    if barcode not in self.full_expression_profile[gene]:
                        self.full_expression_profile[gene][barcode] = self.exon_expression_profile[gene][barcode]
                    else:
                        self.full_expression_profile[gene][barcode] = self.full_expression_profile[gene][barcode].union(self.exon_expression_profile[gene][barcode])
                    #self.full_expression_profile[gene][barcode] = self.exon_expression_profile[gene][barcode]
                elif self.has_specific_gene_cell_expression_with_feature(gene,barcode,'intron'):
                    if barcode not in self.full_expression_profile[gene]:
                        self.full_expression_profile[gene][barcode] = self.intron_expression_profile[gene][barcode]
                    else:
                        self.full_expression_profile[gene][barcode] = self.full_expression_profile[gene][barcode].union(self.intron_expression_profile[gene][barcode])
                    #self.full_expression_profile[gene][barcode] = self.intron_expression_profile[gene][barcode]

    def get_all_cells(self):
        return self.all_cell_barcode

    def get_all_genes(self):
        return self.all_genes

    def get_cells_with_feature(self,feature):
        if feature == 'exon':
            return self.exon_cell_barcode
        elif feature == 'intron':
            return self.intron_cell_barcode
        elif feature == 'total':
            return self.all_cell_barcode
        elif feature == 'direct':
            return self.direct_cell_barcode

    def get_genes_with_feature(self,feature):
        if feature == 'exon':
            return self.exon_genes
        elif feature == 'intron':
            return self.intron_genes
        elif feature == 'total':
            return self.all_genes
        elif feature == 'direct':
            return self.direct_genes

    def has_specific_gene_cell_expression_with_feature(self,gene,cell,feature):
        if feature == 'exon' and gene in self.exon_genes:
            return cell in self.exon_expression_profile[gene]
        elif feature == 'intron' and gene in self.intron_genes:
            return cell in self.intron_expression_profile[gene]
        elif feature == 'total' and gene in self.all_genes:
            return cell in self.full_expression_profile[gene]
        elif feature == 'direct' and gene in self.direct_genes:
            return cell in self.direct_deter_full_expression_profile[gene]
        else:
            False

    def get_gene_cell_expression_with_feature(self,gene,cell,feature):
        if feature == 'exon':
            return self.exon_expression_profile[gene][cell]
        elif feature == 'intron':
            return self.intron_expression_profile[gene][cell]
        elif feature == 'total':
            return self.full_expression_profile[gene][cell]
        elif feature == 'direct':
            return self.direct_deter_full_expression_profile[gene][cell]

    def _umi_to_counts_ddict(self,ddict):
        for gene in ddict.keys():
            for cell in ddict[gene].keys():
                ddict[gene][cell]=len(ddict[gene][cell])
        return ddict
