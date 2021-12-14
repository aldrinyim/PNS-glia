#!/usr/bin/python3

'''
Private Class, do not modify
'''
from collections import defaultdict
import numpy as np

class Transcriptional_modules:

    def __init__(self):
        self.entrez_ensembl_dict = {}
        self.entrez_pathway_dict = {}
        self.pathway_entrez_dict = {}
        self.map_name_dict = {}
        #self.all_maps = []

    def get_all_modules(self):
        return list(self.pathway_entrez_dict.keys())
        #print(self.all_maps)
        #return self.all_maps

    def put_entrez_ensembl(self,address):
        with open(address) as self.f:
            for self.line in self.f:
                self.data=self.line.rstrip().split('\t')
                self.entrez = self.data[0].split(':')[1]
                self.ensembl = self.data[1].split(':')[1]
                self.entrez_ensembl_dict[self.entrez] = self.ensembl
        self.f.close()

    def put_entrez_pathway(self,address):
        #Check gene if exist in entrez database
        with open(address) as self.f:
            for self.line in self.f:
                self.data = self.line.rstrip().split('\t')
                self.temp_entrez = self.data[0].split(':')[1]
                self.temp_map = self.data[1][8:]

                if self.temp_entrez in self.entrez_ensembl_dict.keys():
                    if self.temp_entrez not in self.entrez_pathway_dict:
                        self.entrez_pathway_dict[self.temp_entrez] = [self.temp_map]
                    else:
                        self.entrez_pathway_dict[self.temp_entrez].append(self.temp_map)

                    if self.temp_map not in self.pathway_entrez_dict:
                        self.pathway_entrez_dict[self.temp_map]=[self.temp_entrez]
                    else:
                        self.pathway_entrez_dict[self.temp_map].append(self.temp_entrez)
        self.f.close()

    def put_map_name(self,address):
        with open(address) as self.f:
            for self.line in self.f:
                self.data = self.line.rstrip().split('\t')
                self.map_name_dict[self.data[0]]=self.data[1]
        self.f.close()

    def get_map_name(self,id):
        return self.map_name_dict[id]

    def get_gene_by_pathwayID(self,pathwayID):
        return self.pathway_entrez_dict[pathwayID]

    def get_ensembl_by_pathwayID(self,pathwayID):
        return [self.entrez_ensembl_dict[x] for x in self.pathway_entrez_dict[pathwayID] ]
