#!/usr/bin/python3

'''
Private Class, do not modify
'''
from collections import defaultdict
from transcriptional_modules import Transcriptional_modules
from gtf import Gtf
import numpy as np
import time
from tqdm import tqdm

class ScProfiles:

    def __init__(self):
        self.cellExpressionProfile = defaultdict(dict)
        self.cellModulesProfile = defaultdict(lambda: defaultdict(dict))
        self.cellHeader = []
        self.geneName = []
        self.clusterdict = defaultdict(list)
        self.cluster = set()
        #self.tgtf = Gtf()
        #self.ttm = Transcriptional_modules()

    def put_cluster_info(self,paired_list):
        self.clusterdict[int(paired_list[1])].append(paired_list[0])
        #self.clusterdict[paired_list[0]]=paired_list[1]
        self.cluster.add(int(paired_list[1]))

    def get_clusters(self):
        return list(self.cluster)

    def get_exp_cell_gene(self,cell,gene):
        return self.cellExpressionProfile[cell][gene]

    def define_cells(self,data):
        self.cellHeader = data
        for self.name in self.cellHeader:
            self.cellExpressionProfile[self.name] = {}

    def put_gene_for_cells(self,data):
        self.t_gene = data[0]
        for self.i in range(1,len(data)):
            self.t_cellName = self.cellHeader[self.i]
            #put expression to cell and gene
            if self.t_cellName in self.cellExpressionProfile.keys():
                self.cellExpressionProfile[self.t_cellName][self.t_gene]=float(data[self.i])
            else:
                print("Something is wrong")

    def convert_ensembl_to_genename(tlist):
        return [mm10_gtf.get_gene_name(x) for x in tlist]

    def define_cell_modules(self,ttm,tgtf):
        #self.tgtf = tgtf
        #self.ttm = ttm
        self.moduleIDs = ttm.get_all_modules()
        print("Number of modules : ",len(self.moduleIDs))
        print("Total number of cells : ", len(self.cellHeader)-1)
        self.lost_counter = 0
        self.cell_counter = 0
        self.module_counter = 0
        self.curr_module = ""

        for self.i in range(1,len(self.cellHeader)):
            self.name = self.cellHeader[self.i]
            self.cell_counter+=1
            if self.cell_counter % 100 == 0:
                print("Number of cells processed : ", self.cell_counter, end="\r")

            for self.module in self.moduleIDs:
                self.ensembls = ttm.get_ensembl_by_pathwayID(self.module)
                self.genes = [tgtf.get_gene_name(self.x) for self.x in self.ensembls]
                for self.gene in self.genes:
                    if self.gene is not None:
                        if self.gene in self.cellExpressionProfile[self.name].keys():
                            self.cellModulesProfile[self.name][self.module][self.gene] = self.cellExpressionProfile[self.name][self.gene]
                        else:
                            self.cellModulesProfile[self.name][self.module][self.gene] = 0.0
                    else:
                        self.lost_counter+=1
        print("")

    #def get_cell_module_expression_per_gene(self,)

    def get_expression_modules_by_cluster(self,clusterID,moduleID):
        self.cells_in_cluster = self.clusterdict[clusterID]
        self.group_expression_modules = []
        for self.cell in self.cells_in_cluster:
            self.group_expression_modules.append(self.cellModulesProfile[self.cell][moduleID])
        return self.group_expression_modules

    #def get_cellName_by_cluster(self,clusterID):
    #    return filter(lambda x: x in self.clusterdict[clusterID], )

    def has_gene(self,cell,geneName):
        if geneName in self.cellExpressionProfile[cell].keys():
            return True
        else:
            return False

    def get_cell_expression_profile_by_moduleID(self,cell,moduleID):
        return self.cellModulesProfile[cell][moduleID]
