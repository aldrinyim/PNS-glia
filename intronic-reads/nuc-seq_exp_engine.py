#!/usr/bin/python3
"""
Author: Aldrin Yim
Version 0.2 (beta)
Date: Feb 2018
Usage: python3 nuc-seq_exp_engine.py --help

Note1: Support Multiprocessing!
"""
import sys, os, subprocess, argparse, pysam, pandas as pd, seaborn as sns, matplotlib.pyplot as plt, math, statistics, numpy as np
from collections import ChainMap, defaultdict
import multiprocessing, tqdm
from gtf import Gtf
from transcriptome_profile import Transcriptome_profile
from operator import itemgetter
from scipy import stats

nuclei_profile = Transcriptome_profile()

current_wd = os.getcwd()
print('working directory: ',current_wd)

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--bam", help=" ESSENTIAL Input BAM file")
parser.add_argument("-g", "--gtf", help=" ESSENTIAL GTF file")
parser.add_argument("-c", "--core", help=" Support multi-core")
parser.add_argument("-a", "--analysis", help=" Analysis, model fitting and correction")
parser.add_argument("-s", "--selected", help=" Start analysis with selected genes")
parser.add_argument("-e", "--express", help=" Pick top N expressed genes for modeling")
parser.add_argument("-n", "--no", help=" When run on HTCF, disable graph drawing with seaborn and matplotlib")
parser.add_argument("-d", "--debug", help=" Testing experimental functions, debug only")
parser.add_argument("-y", "--coverageAnalysis", help=" Pick top N expressed genes to perform Coverage analysis on intron and exon regions separately")
parser.add_argument("-w", "--writeMatrix", help=" Write single nuclei expression matrix (10X format) for downstream analysis")
parser.add_argument("-r", "--directGeneOnly", help="Determine unique umi count for full gene body, ignore exon/intron count")
parser.add_argument("-m", "--multiplebam", help="Ad-hoc, only for 3 bam files at one time")
args = parser.parse_args()

if args.bam == None or args.gtf == None:
    parser.print_help()
    sys.exit()

if args.core == None:
    args.core = 1
else:
    args.core = int(args.core)

def remove_quote(word):
    return word.split('\'')[1]

def full_gene_unique_umi_worker(target_gene):
    expression_dict = defaultdict(dict)
    bamfile = pysam.AlignmentFile(args.bam,"rb")
    i_info = mm10_gtf.get_gene_info(target_gene)
    #print(i_info[2],'|',i_info[3],'|',i_info[4])
    read_stack = bamfile.fetch(i_info[2],i_info[3],i_info[4])
    for read in read_stack:
        #print(read.mapping_quality)
        if read.mapping_quality == 255:
            #print(read)
            #print(read.mapping_quality)
            read_tags = list(read.tags)
            map_type = ''
            cellular_id = ''
            umi = ''
            for indiv_tag in read_tags:
                tag_name = indiv_tag[0]
                tag_value = indiv_tag[1]
                if tag_name == 'CB':
                    cellular_id = tag_value
                if tag_name == 'UB':
                    umi = tag_value
                if tag_name == 'RE':
                    map_type = tag_value
            #print('RE tag = ', map_type)
            if cellular_id != '' and umi != '' and len(umi) >= 9:
                if cellular_id not in expression_dict[target_gene]:
                    expression_dict[target_gene][cellular_id]=set()
                    expression_dict[target_gene][cellular_id].add(umi)
                else:
                    expression_dict[target_gene][cellular_id].add(umi)
    #print(expression_dict)
    return _umi_to_counts_ddict(expression_dict)

def test_full_gene_unique_umi_worker_multiple_bam_unlimited(target_gene):
    expression_dict = defaultdict(dict)
    i_info = mm10_gtf.get_gene_info(target_gene)
    bamfiles = str(args.bam).split(',')

    bamstream_list = []
    for i in bamfiles:
        bamstream_list.append(pysam.AlignmentFile(i,"rb"))

    readstream_list = []
    for j in bamstream_list:
        readstream_list.append(j.fetch(i_info[2],i_info[3],i_info[4]))

    for k in range(0,len(readstream_list)):
        for read in readstream_list[k]:
            if read.mapping_quality == 255:
                #print(read)
                #print(read.mapping_quality)
                read_tags = list(read.tags)
                map_type = ''
                cellular_id = ''
                umi = ''
                for indiv_tag in read_tags:
                    tag_name = indiv_tag[0]
                    tag_value = indiv_tag[1]
                    if tag_name == 'CB':
                        cellular_id = tag_value
                    if tag_name == 'UB':
                        umi = tag_value
                    if tag_name == 'RE':
                        map_type = tag_value
                #print('RE tag = ', map_type)
                if cellular_id != '' and umi != '' and len(umi) >= 9:
                    cellular_id = cellular_id+str(k)
                    if cellular_id not in expression_dict[target_gene]:
                        expression_dict[target_gene][cellular_id]=set()
                        expression_dict[target_gene][cellular_id].add(umi)
                    else:
                        expression_dict[target_gene][cellular_id].add(umi)

    for j in bamstream_list:
        j.close()

    return _umi_to_counts_ddict(expression_dict)

def full_gene_unique_umi_worker_multiple_bam(target_gene):
    expression_dict = defaultdict(dict)
    bamfiles = str(args.bam).split(',')
    #print('bamfiles input : ',bamfiles)
    bamfile1 = pysam.AlignmentFile(bamfiles[0],"rb")
    bamfile2 = pysam.AlignmentFile(bamfiles[1],"rb")
    bamfile3 = pysam.AlignmentFile(bamfiles[2],"rb")
    bamfile4 = pysam.AlignmentFile(bamfiles[3],"rb")
    i_info = mm10_gtf.get_gene_info(target_gene)
    #print(i_info[2],'|',i_info[3],'|',i_info[4])
    read_stack1 = bamfile1.fetch(i_info[2],i_info[3],i_info[4])
    read_stack2 = bamfile2.fetch(i_info[2],i_info[3],i_info[4])
    read_stack3 = bamfile3.fetch(i_info[2],i_info[3],i_info[4])
    read_stack4 = bamfile4.fetch(i_info[2],i_info[3],i_info[4])
    for read in read_stack1:
        #print(read.mapping_quality)
        if read.mapping_quality == 255:
            #print(read)
            #print(read.mapping_quality)
            read_tags = list(read.tags)
            map_type = ''
            cellular_id = ''
            umi = ''
            for indiv_tag in read_tags:
                tag_name = indiv_tag[0]
                tag_value = indiv_tag[1]
                if tag_name == 'CB':
                    cellular_id = tag_value
                if tag_name == 'UB':
                    umi = tag_value
                if tag_name == 'RE':
                    map_type = tag_value
            #print('RE tag = ', map_type)
            if cellular_id != '' and umi != '' and len(umi) >= 9:
                cellular_id = cellular_id+'1'
                if cellular_id not in expression_dict[target_gene]:
                    expression_dict[target_gene][cellular_id]=set()
                    expression_dict[target_gene][cellular_id].add(umi)
                else:
                    expression_dict[target_gene][cellular_id].add(umi)
    for read in read_stack2:
        #print(read.mapping_quality)
        if read.mapping_quality == 255:
            #print(read)
            #print(read.mapping_quality)
            read_tags = list(read.tags)
            map_type = ''
            cellular_id = ''
            umi = ''
            for indiv_tag in read_tags:
                tag_name = indiv_tag[0]
                tag_value = indiv_tag[1]
                if tag_name == 'CB':
                    cellular_id = tag_value
                if tag_name == 'UB':
                    umi = tag_value
                if tag_name == 'RE':
                    map_type = tag_value
            #print('RE tag = ', map_type)
            if cellular_id != '' and umi != '' and len(umi) >= 9:
                cellular_id = cellular_id+'2'
                if cellular_id not in expression_dict[target_gene]:
                    expression_dict[target_gene][cellular_id]=set()
                    expression_dict[target_gene][cellular_id].add(umi)
                else:
                    expression_dict[target_gene][cellular_id].add(umi)
    for read in read_stack3:
        #print(read.mapping_quality)
        if read.mapping_quality == 255:
            #print(read)
            #print(read.mapping_quality)
            read_tags = list(read.tags)
            map_type = ''
            cellular_id = ''
            umi = ''
            for indiv_tag in read_tags:
                tag_name = indiv_tag[0]
                tag_value = indiv_tag[1]
                if tag_name == 'CB':
                    cellular_id = tag_value
                if tag_name == 'UB':
                    umi = tag_value
                if tag_name == 'RE':
                    map_type = tag_value
            #print('RE tag = ', map_type)
            if cellular_id != '' and umi != '' and len(umi) >= 9:
                cellular_id = cellular_id+'3'
                if cellular_id not in expression_dict[target_gene]:
                    expression_dict[target_gene][cellular_id]=set()
                    expression_dict[target_gene][cellular_id].add(umi)
                else:
                    expression_dict[target_gene][cellular_id].add(umi)
    for read in read_stack4:
        #print(read.mapping_quality)
        if read.mapping_quality == 255:
            #print(read)
            #print(read.mapping_quality)
            read_tags = list(read.tags)
            map_type = ''
            cellular_id = ''
            umi = ''
            for indiv_tag in read_tags:
                tag_name = indiv_tag[0]
                tag_value = indiv_tag[1]
                if tag_name == 'CB':
                    cellular_id = tag_value
                if tag_name == 'UB':
                    umi = tag_value
                if tag_name == 'RE':
                    map_type = tag_value
            #print('RE tag = ', map_type)
            if cellular_id != '' and umi != '' and len(umi) >= 9:
                cellular_id = cellular_id+'4'
                if cellular_id not in expression_dict[target_gene]:
                    expression_dict[target_gene][cellular_id]=set()
                    expression_dict[target_gene][cellular_id].add(umi)
                else:
                    expression_dict[target_gene][cellular_id].add(umi)
    #print(expression_dict)
    bamfile1.close()
    bamfile2.close()
    bamfile3.close()
    bamfile4.close()
    return _umi_to_counts_ddict(expression_dict)

def exon_worker(target_gene):
    #print('In worker, ',target_gene)
    #temp_nuclei_exp = Transcriptome_profile()
    bamfile = pysam.AlignmentFile(args.bam,"rb")
    i_info = mm10_gtf.get_gene_info(target_gene)

    exon_info = mm10_gtf.get_exon_map_specific_gene(target_gene)
    #intron_info = mm10_gtf.get_intron_map_specific_gene(target_gene)
    #print('exon info - ',exon_info)
    #print('intron info - ', intron_info)
    gene_chrom = mm10_gtf.get_gene_chrom(target_gene)
    #eexon_exp_level = []
    eexon_ind_exp = []
    gene_based_exon_count = set()
    #if target_gene == "ENSMUSG00000025900":
    #print(target_gene,'|',i_info[2],'|',exon_info)
    for eexon in exon_info:
        #print(target_gene,' exon |',eexon)
        #temp_exon_umi = set()
        temp_reads = 0
        read_stack = bamfile.fetch(i_info[2],eexon[0],eexon[1])
        for read in read_stack:
            if read.mapping_quality == 255:
                #print(read)
                #print(read.mapping_quality)
                read_tags = list(read.tags)
                map_type = ''
                cellular_id = ''
                umi = ''
                for indiv_tag in read_tags:
                    tag_name = indiv_tag[0]
                    tag_value = indiv_tag[1]
                    if tag_name == 'CB':
                        cellular_id = tag_value
                    if tag_name == 'UB':
                        umi = tag_value
                    if tag_name == 'RE':
                        map_type = tag_value
                    if tag_name != '' and cellular_id != '' and len(umi) >= 9:
                        #temp_exon_umi.add(umi)
                        gene_based_exon_count.add(umi)
                        temp_reads+=1
            #eexon_exp_level.append(len(temp_exon_umi))
        eexon_ind_exp.append([gene_chrom,eexon[0],eexon[1],temp_reads])
    #if target_gene == "ENSMUSG00000025900":
    #    print('ENSMUSG00000025900 : ', gene_based_exon_count)
    bamfile.close()
    #return [{target_gene:sum(eexon_exp_level)},{target_gene:len(gene_based_exon_count)},{target_gene:eexon_ind_exp}]
    return [{target_gene:len(gene_based_exon_count)},{target_gene:eexon_ind_exp}]

def intron_worker(target_gene):
    #print('In worker, ',target_gene)
    #temp_nuclei_exp = Transcriptome_profile()
    bamfile = pysam.AlignmentFile(args.bam,"rb")
    i_info = mm10_gtf.get_gene_info(target_gene)

    #exon_info = mm10_gtf.get_exon_map_specific_gene(target_gene)
    intron_info = mm10_gtf.get_intron_map_specific_gene(target_gene)
    gene_chrom = mm10_gtf.get_gene_chrom(target_gene)
    #intron_exp_level = []
    iintron_ind_exp = []
    gene_based_intron_count = set()
    for eintron in intron_info:
        #print(target_gene,' intron |',eintron)
        #temp_intron_umi = set()
        temp_reads = 0
        read_stack = bamfile.fetch(i_info[2],eintron[0],eintron[1])
        for read in read_stack:
            if read.mapping_quality == 255:
                #print(read)
                #print(read.mapping_quality)
                read_tags = list(read.tags)
                map_type = ''
                cellular_id = ''
                umi = ''
                for indiv_tag in read_tags:
                    tag_name = indiv_tag[0]
                    tag_value = indiv_tag[1]
                    if tag_name == 'CB':
                        cellular_id = tag_value
                    if tag_name == 'UB':
                        umi = tag_value
                    if tag_name == 'RE':
                        map_type = tag_value
                    if tag_name != '' and cellular_id != '' and len(umi) >= 9:
                        #temp_intron_umi.add(umi)
                        gene_based_intron_count.add(umi)
                        temp_reads+=1
            #intron_exp_level.append(len(temp_intron_umi))
        iintron_ind_exp.append([gene_chrom,eintron[0],eintron[1],temp_reads])
    #print(target_gene:[sum(eexon_exp_level),sum(intron_exp_level)])
    bamfile.close()
    #return [{target_gene:sum(intron_exp_level)},{target_gene:len(gene_based_intron_count)}]
    return [{target_gene:len(gene_based_intron_count)},{target_gene:iintron_ind_exp}]

def single_nuclei_exon_intron_worker(target_gene,features):
    expression_dict = defaultdict(dict)
    bamfile = pysam.AlignmentFile(args.bam,"rb")
    i_info = mm10_gtf.get_gene_info(target_gene)

    if features == "intron":
        info = mm10_gtf.get_intron_map_specific_gene(target_gene)
    elif features == "exon":
        info = mm10_gtf.get_exon_map_specific_gene(target_gene)

    for pos in info:
        read_stack = bamfile.fetch(i_info[2],pos[0],pos[1])
        for read in read_stack:
            #print(read.mapping_quality)
            if read.mapping_quality == 255:
                #print(read)
                #print(read.mapping_quality)
                read_tags = list(read.tags)
                map_type = ''
                cellular_id = ''
                umi = ''
                for indiv_tag in read_tags:
                    tag_name = indiv_tag[0]
                    tag_value = indiv_tag[1]
                    if tag_name == 'CB':
                        cellular_id = tag_value
                    if tag_name == 'UB':
                        umi = tag_value
                    if tag_name == 'RE':
                        map_type = tag_value
                #print('RE tag = ', map_type)
                if cellular_id != '' and umi != '' and len(umi) >= 9:
                    if cellular_id not in expression_dict[target_gene]:
                        expression_dict[target_gene][cellular_id]=set()
                        expression_dict[target_gene][cellular_id].add(umi)
                    else:
                        expression_dict[target_gene][cellular_id].add(umi)
    return expression_dict
    #return _umi_to_counts_ddict(expression_dict)

def _umi_to_counts_ddict(ddict):
    for gene in ddict.keys():
        for cell in ddict[gene].keys():
            ddict[gene][cell]=len(ddict[gene][cell])
    return ddict

'''
def intron_pileup_worker(target_gene):
    bamfile = pysam.AlignmentFile(args.bam,"rb")
    i_info = mm10_gtf.get_gene_info(target_gene)

    gene_chrom = mm10_gtf.get_gene_chrom(target_gene)
    gene_start = mm10_gtf.get_gene_start(target_gene)
    gene_end = mm10_gtf.get_gene_end(target_gene)
    gene_length = gene_end-gene_start+1

    exon_info = mm10_gtf.get_exon_map_specific_gene(target_gene)
    intron_info = mm10_gtf.get_intron_map_specific_gene(target_gene)

    exon_counts = []
    intron_counts = []

    for each_block in exon_info:
        pileup_stack = bamfile.pileup(gene_chrom,each_block[0],each_block[1]+1,max_depth=100000000,truncate=True)
'''

def define_graphic_region(target_gene):
    #print("In define_graphic_region, gene : ",target_gene)
    bamfile = pysam.AlignmentFile(args.bam,"rb")
    i_info = mm10_gtf.get_gene_info(target_gene)

    gene_chrom = mm10_gtf.get_gene_chrom(target_gene)
    gene_start = mm10_gtf.get_gene_start(target_gene)
    gene_end = mm10_gtf.get_gene_end(target_gene)
    gene_length = gene_end-gene_start+1
    gene_strand = mm10_gtf.get_gene_strand(target_gene)
    #print('gene strand : ',gene_strand)
    exon_info = mm10_gtf.get_exon_map_specific_gene(target_gene)
    intron_info = mm10_gtf.get_intron_map_specific_gene(target_gene)
    #print('Gene start | End - ', gene_start,'|',gene_end)
    #print('Exon map : ', exon_info)
    #print('Intron map : ', intron_info)

    exon_length = _get_length(exon_info)
    intron_length = _get_length(intron_info)

    if exon_length >= 301 and intron_length >= 301:
        exon_cov = []
        intron_cov = []

        ##Construct full gene body coverage first, Current version does not consider UMI
        #gene_cov = [0]*gene_length
        start_pos = gene_start
        pileup_stack = bamfile.pileup(gene_chrom,gene_start,gene_end+1,max_depth=100000000,truncate=True)
        for puc in pileup_stack:

            curr_pos = puc.pos
            curr_cov = puc.n
            exon_checker = _check_pos_in_block(curr_pos,exon_info)
            intron_checker = _check_pos_in_block(curr_pos,intron_info)

            while curr_pos != start_pos:        ##Beginning, check for zero coverage
                #print("Coverage = 0, curr pos = ", curr_pos)
                if exon_checker[0] == True:
                    exon_cov.append(0)
                elif intron_checker[0] == True:
                    intron_cov.append(0)
                start_pos+=1

            if exon_checker[0] == True:
                exon_cov.append(curr_cov)
            elif intron_checker[0] == True:
                intron_cov.append(curr_cov)

            start_pos+=1
        bamfile.close()

        if len(exon_cov) == 0 or len(intron_cov) == 0:
            return False
        #print('Exon len by default : ', exon_length)
        #print('Exon len by cov : ', len(exon_cov))
        #print('intron len by default : ', intron_length)
        #print('intron len by cov : ', len(intron_cov))

        max_exon_cov = max(exon_cov)
        max_intron_cov = max(intron_cov)

        if gene_strand == '-':                  #Correct for strand
            exon_cov.reverse()
            intron_cov.reverse()

        #exon_unit_len = _get_percentile_unit(exon_length)
        #intron_unit_len = _get_percentile_unit(intron_length)

        #For exon
        #print('Exon Coverage : ',exon_cov)
        exon_coverage_in_percentile = _get_coverage_in_percentile(exon_cov)
        intron_coverage_in_percentile = _get_coverage_in_percentile(intron_cov)

        return [_normalize_coverage(exon_coverage_in_percentile),_normalize_coverage(intron_coverage_in_percentile)]
    else:
        return False

def _combine_normalized_coverage_percentile(list_of_list_of_coverage):
    return list(map(statistics.mean, zip(*list_of_list_of_coverage)))

def _get_coverage_in_percentile(ori_coverage_list):
    coverage_in_percentile = [0]*100
    unit_length = _get_percentile_unit(len(ori_coverage_list))
    #print(len(ori_coverage_list),'|',unit_length)
    for i in range(0,100):
        if i == 0:
            unit_start = 0
            unit_end = unit_start + unit_length
        elif i > 0:
            unit_start = unit_end + 1
            unit_end = unit_start + unit_length - 1
        elif i == 99:
            unit_start = unit_end + 1
            unit_end = len(ori_coverage_list)
        #print(i,'|',unit_start,'|',unit_end)
        if len(ori_coverage_list[unit_start:unit_end+1]) > 1:
            coverage_in_percentile[i]=statistics.mean(ori_coverage_list[unit_start:unit_end+1])
        else:
            coverage_in_percentile[i]=ori_coverage_list[unit_start:unit_end+1]
    #print("In Function : ")
    #print(coverage_in_percentile)
    #print("End function")
    return coverage_in_percentile

def _normalize_coverage(percentile_coverage_list):
    return list(x/max(percentile_coverage_list) for x in percentile_coverage_list)

def _check_pos_in_block(pos,list_of_coors):
    checker = (False,[])
    for each_block in list_of_coors:
        if pos >= each_block[0] and pos <= each_block[1]:
            checker = (True,each_block)
    return checker

def _get_length(list_of_coors):
    distance = 0
    for each_block in list_of_coors:
        distance += each_block[1]-each_block[0]+1
    return distance

def _get_percentile_unit(length):
    return math.floor(length/100)

def merging_exon_intron_expression_dict(from_exons_worker):         ##Not suggest to use, maybe slower
    gs_dict = from_exons_worker[0]
    for i in range(1,len(from_exons_worker)):
        gs_dict.update(from_exons_worker[i])
    return gs_dict
    #sorted(gs_dict, key=gs.get, reverse=True)

def merging_exon_intron_expression_dict_CM(from_workers):
    return dict(ChainMap(*from_workers))

def merging_single_nuclei_exon_intron_umi_dict_CM(from_workers):
    # pretty crazy, from workers are default dictionary
    # convert to dict for each, iterate, chainmap and dict again
    return dict(ChainMap(*(dict(x) for x in from_workers)))

def merging_multiple_defaultdict(final_profile, list_of_dicts):
    #Assuming no crash in data, as each core is responsible for one gene, and the structure should not crash
    if args.analysis is None:
        total_dicts = len(list_of_dicts)
        ground_state = list_of_dicts[0].get_profile()
        for i in range(1,len(list_of_dicts)):
            ground_state.update(list_of_dicts[i].get_profile())
        final_profile.reset_profile(ground_state)
    else:
        total_dicts = len(list_of_dicts)
        ground_state = list_of_dicts[0].get_profile()
        exon_ground_state = list_of_dicts[0].get_exon_profile()
        intron_ground_state = list_of_dicts[0].get_intron_profile()
        for i in range(1,len(list_of_dicts)):
            ground_state.update(list_of_dicts[i].get_profile())
            exon_ground_state.update(list_of_dicts[i].get_exon_profile())
            intron_ground_state.update(list_of_dicts[i].get_intron_profile())
        final_profile.reset_profile(ground_state)
        final_profile.reset_exon_profile(exon_ground_state)
        final_profile.reset_intron_profile(intron_ground_state)


mm10_gtf = Gtf()
with open(args.gtf) as f:
    for line in f:
        readin=line.rstrip()
        if readin[0] != '#':
            mm10_gtf.put_gene(readin)
f.close()

total_genes = mm10_gtf.get_total_ids()
print('Curating all exon info, no of genes - ',len(total_genes))

mm10_gtf.create_intron_map()
mm10_gtf.create_exon_map()

print("HHV1gp067" in total_genes)

#intron_map = mm10_gtf.get_intron_map()
#exon_map = mm10_gtf.get_exon_map()

manager = multiprocessing.Manager()
lock = manager.Lock()

if args.directGeneOnly is None and args.multiplebam is None:
    with multiprocessing.Pool(processes=args.core) as p:
        exon_outputs = list(tqdm.tqdm(p.imap(exon_worker,total_genes),total=len(total_genes)))
    p.close()
    p.join()

    d=merging_exon_intron_expression_dict_CM(list(list(zip(*exon_outputs))[0]))  ##Using unique exon count
    #print('After merging ; ',d['ENSMUSG00000025900'])
    outf = open('gene_expression_by_exon.txt','w')
    sorted_gene_expression_by_exon = [(k, d[k]) for k in sorted(d, key=d.get, reverse=True)]
    for k, v in sorted_gene_expression_by_exon:
        outf.write(k+'\t'+mm10_gtf.get_gene_name(k)+'\t'+str(v)+'\n')
    outf.close()

    exd = merging_exon_intron_expression_dict_CM(list(list(zip(*exon_outputs))[1]))
    top_5000_genes = sorted_gene_expression_by_exon[0:5000]
    outf = open('exon_based_individual_expression.txt','w')
    outf.write('#Zero-based_coordinate_system\n')
    outf.write('gene_id\tgene_name\tstrand\tchromosome\texon_start\texon_end\texon_unique_read_count\n')
    for i,v in top_5000_genes:
        #outf.write(i+'\t'+mm10_gtf.get_gene_name(i)+'\t'+mm10_gtf.get_gene_strand(i)+'\t')
        exon_blocks = exd[i]
        for block in exon_blocks:
            outf.write(i+'\t'+mm10_gtf.get_gene_name(i)+'\t'+mm10_gtf.get_gene_strand(i)+'\t'+str(block[0])+'\t'+str(block[1])+'\t'+str(block[2])+'\t'+str(block[3])+'\n')
            #outf.write(i+'\t'+mm10_gtf.get_gene_name(i)+'\t'+str(exd[i][0])+'\t'+str(exd[i][1])+'\t'+str(exd[i][2])+'\t'+str(mm10_gtf.get_gene_strand(i))+'\t'+str(exd[i][0])+'\n')
    outf.close()

if args.directGeneOnly is not None and args.multiplebam is not None:
    #bamfiles = str(args.bam).split(',')
    #print(bamfiles)
    with multiprocessing.Pool(processes=args.core) as p:
        #full_gene_outputs = list(tqdm.tqdm(p.imap(full_gene_unique_umi_worker_multiple_bam,total_genes),total=len(total_genes)))
        full_gene_outputs = list(tqdm.tqdm(p.imap(test_full_gene_unique_umi_worker_multiple_bam_unlimited,total_genes),total=len(total_genes)))
    p.close()
    p.join()

    full_gene_dict = merging_single_nuclei_exon_intron_umi_dict_CM(full_gene_outputs)
    nuclei_profile.put_full_gene_profile(full_gene_dict)
    #print(_umi_to_counts_ddict(full_gene_dict.copy()))

    direct_cells = list(nuclei_profile.get_cells_with_feature('direct'))
    direct_genes = list(nuclei_profile.get_genes_with_feature('direct'))

    if not os.path.exists(current_wd+'/direct_expression/'):
        os.makedirs(current_wd+'/direct_expression/')

    barcodes_outf = open(current_wd+'/direct_expression/barcodes.tsv','w')
    for cell in direct_cells:
        barcodes_outf.write(cell+'\n')
    barcodes_outf.close()

    genes_outf = open(current_wd+'/direct_expression/genes.tsv','w')
    for gene in direct_genes:
        genes_outf.write(gene+'\t'+mm10_gtf.get_gene_name(gene)+'\n')
    genes_outf.close()

    combinations = 0
    for i in range(0,len(direct_genes)):
        for j in range(0,len(direct_cells)):
            if nuclei_profile.has_specific_gene_cell_expression_with_feature(direct_genes[i],direct_cells[j],'direct'):
                combinations+=1

    matrix_outf = open(current_wd+'/direct_expression/matrix.mtx','w')
    matrix_outf.write('%%MatrixMarket matrix coordinate integer general\n%\n')
    matrix_outf.write(str(len(direct_genes))+' '+str(len(direct_cells))+' '+str(combinations)+'\n')
    for i in range(0,len(direct_genes)):
        for j in range(0,len(direct_cells)):
            if nuclei_profile.has_specific_gene_cell_expression_with_feature(direct_genes[i],direct_cells[j],'direct'):
                matrix_outf.write(str(i+1)+' '+str(j+1)+' '+str(nuclei_profile.get_gene_cell_expression_with_feature(direct_genes[i],direct_cells[j],'direct'))+'\n')
    matrix_outf.close()
    #print('End of exon profiling, ', exon_outputs)

if args.directGeneOnly is not None and args.multiplebam is None:
    with multiprocessing.Pool(processes=args.core) as p:
        full_gene_outputs = list(tqdm.tqdm(p.imap(full_gene_unique_umi_worker,total_genes),total=len(total_genes)))
    p.close()
    p.join()

    full_gene_dict = merging_single_nuclei_exon_intron_umi_dict_CM(full_gene_outputs)
    nuclei_profile.put_full_gene_profile(full_gene_dict)
    #print(_umi_to_counts_ddict(full_gene_dict.copy()))

    direct_cells = list(nuclei_profile.get_cells_with_feature('direct'))
    direct_genes = list(nuclei_profile.get_genes_with_feature('direct'))

    if not os.path.exists(current_wd+'/direct_expression/'):
        os.makedirs(current_wd+'/direct_expression/')

    barcodes_outf = open(current_wd+'/direct_expression/barcodes.tsv','w')
    for cell in direct_cells:
        barcodes_outf.write(cell+'\n')
    barcodes_outf.close()

    genes_outf = open(current_wd+'/direct_expression/genes.tsv','w')
    for gene in direct_genes:
        genes_outf.write(gene+'\t'+mm10_gtf.get_gene_name(gene)+'\n')
    genes_outf.close()

    combinations = 0
    for i in range(0,len(direct_genes)):
        for j in range(0,len(direct_cells)):
            if nuclei_profile.has_specific_gene_cell_expression_with_feature(direct_genes[i],direct_cells[j],'direct'):
                combinations+=1

    matrix_outf = open(current_wd+'/direct_expression/matrix.mtx','w')
    matrix_outf.write('%%MatrixMarket matrix coordinate integer general\n%\n')
    matrix_outf.write(str(len(direct_genes))+' '+str(len(direct_cells))+' '+str(combinations)+'\n')
    for i in range(0,len(direct_genes)):
        for j in range(0,len(direct_cells)):
            if nuclei_profile.has_specific_gene_cell_expression_with_feature(direct_genes[i],direct_cells[j],'direct'):
                matrix_outf.write(str(i+1)+' '+str(j+1)+' '+str(nuclei_profile.get_gene_cell_expression_with_feature(direct_genes[i],direct_cells[j],'direct'))+'\n')
    matrix_outf.close()
    #print('End of exon profiling, ', exon_outputs)

if args.writeMatrix is not None:

    print('About to write Matrix, profile all introns')

    with multiprocessing.Pool(processes=args.core) as p:
        exon_single_nuc_outputs = list(tqdm.tqdm(p.starmap(single_nuclei_exon_intron_worker,list(list(x) for x in list(zip(total_genes,['exon']*len(total_genes))))),total=len(total_genes)))
        #exon_single_nuc_outputs = list(tqdm.tqdm(p.map(single_nuclei_exon_intron_worker,total_genes,['exon']*len(total_genes)),total=len(total_genes)))
    p.close()
    p.join()


    with multiprocessing.Pool(processes=args.core) as p:
        #intron_single_nuc_outputs = list(tqdm.tqdm(p.map(single_nuclei_exon_intron_worker,list(zip(total_genes,['intron']*len(total_genes))),total=len(total_genes))))
        intron_single_nuc_outputs = list(tqdm.tqdm(p.starmap(single_nuclei_exon_intron_worker,list(list(x) for x in list(zip(total_genes,['intron']*len(total_genes))))),total=len(total_genes)))
    p.close()
    p.join()

    print('End of intron profiling ')

    full_exon_single_nuc_outputs = merging_single_nuclei_exon_intron_umi_dict_CM(exon_single_nuc_outputs)
    full_intron_single_nuc_output = merging_single_nuclei_exon_intron_umi_dict_CM(intron_single_nuc_outputs)

    #print('exon data : ',full_exon_single_nuc_outputs)
    #print('intron data : ',full_intron_single_nuc_output)

    nuclei_profile.put_exon_expression_profile(full_exon_single_nuc_outputs)
    nuclei_profile.put_intron_expresion_profile(full_intron_single_nuc_output)
    nuclei_profile.combine_expression_profile()

    all_cells = list(nuclei_profile.get_all_cells())
    all_genes = list(nuclei_profile.get_all_genes())

    #Write MEX format table
    if not os.path.exists(current_wd+'/combined_expression/'):
        os.makedirs(current_wd+'/combined_expression/')

    barcodes_outf = open(current_wd+'/combined_expression/barcodes.tsv','w')
    for cell in all_cells:
        barcodes_outf.write(cell+'\n')
    barcodes_outf.close()

    genes_outf = open(current_wd+'/combined_expression/genes.tsv','w')
    for gene in all_genes:
        genes_outf.write(gene+'\t'+mm10_gtf.get_gene_name(gene)+'\n')
    genes_outf.close()

    combinations = 0
    for i in range(0,len(all_genes)):
        for j in range(0,len(all_cells)):
            if nuclei_profile.has_specific_gene_cell_expression_with_feature(all_genes[i],all_cells[j],'total'):
                combinations+=1

    matrix_outf = open(current_wd+'/combined_expression/matrix.mtx','w')
    matrix_outf.write('%%MatrixMarket matrix coordinate integer general\n%\n')
    matrix_outf.write(str(len(all_genes))+' '+str(len(all_cells))+' '+str(combinations)+'\n')
    for i in range(0,len(all_genes)):
        for j in range(0,len(all_cells)):
            if nuclei_profile.has_specific_gene_cell_expression_with_feature(all_genes[i],all_cells[j],'total'):
                matrix_outf.write(str(i+1)+' '+str(j+1)+' '+str(len(nuclei_profile.get_gene_cell_expression_with_feature(all_genes[i],all_cells[j],'total')))+'\n')
    matrix_outf.close()

    ##For Exon ONLY!
    exon_cells = list(nuclei_profile.get_cells_with_feature('exon'))
    exon_genes = list(nuclei_profile.get_genes_with_feature('exon'))
    if not os.path.exists(current_wd+'/exon_expression/'):
        os.makedirs(current_wd+'/exon_expression/')

    barcodes_outf = open(current_wd+'/exon_expression/barcodes.tsv','w')
    for cell in exon_cells:
        barcodes_outf.write(cell+'\n')
    barcodes_outf.close()

    genes_outf = open(current_wd+'/exon_expression/genes.tsv','w')
    for gene in exon_genes:
        genes_outf.write(gene+'\t'+mm10_gtf.get_gene_name(gene)+'\n')
    genes_outf.close()

    combinations = 0
    for i in range(0,len(exon_genes)):
        for j in range(0,len(exon_cells)):
            if nuclei_profile.has_specific_gene_cell_expression_with_feature(exon_genes[i],exon_cells[j],'exon'):
                combinations+=1

    matrix_outf = open(current_wd+'/exon_expression/matrix.mtx','w')
    matrix_outf.write('%%MatrixMarket matrix coordinate integer general\n%\n')
    matrix_outf.write(str(len(exon_genes))+' '+str(len(exon_cells))+' '+str(combinations)+'\n')
    for i in range(0,len(exon_genes)):
        for j in range(0,len(exon_cells)):
            if nuclei_profile.has_specific_gene_cell_expression_with_feature(exon_genes[i],exon_cells[j],'exon'):
                matrix_outf.write(str(i+1)+' '+str(j+1)+' '+str(len(nuclei_profile.get_gene_cell_expression_with_feature(exon_genes[i],exon_cells[j],'exon')))+'\n')
    matrix_outf.close()

    ##For Introns ONLY!
    intron_cells = list(nuclei_profile.get_cells_with_feature('intron'))
    intron_genes = list(nuclei_profile.get_genes_with_feature('intron'))
    if not os.path.exists(current_wd+'/intron_expression/'):
        os.makedirs(current_wd+'/intron_expression/')

    barcodes_outf = open(current_wd+'/intron_expression/barcodes.tsv','w')
    for cell in intron_cells:
        barcodes_outf.write(cell+'\n')
    barcodes_outf.close()

    genes_outf = open(current_wd+'/intron_expression/genes.tsv','w')
    for gene in intron_genes:
        genes_outf.write(gene+'\t'+mm10_gtf.get_gene_name(gene)+'\n')
    genes_outf.close()

    combinations = 0
    for i in range(0,len(intron_genes)):
        for j in range(0,len(intron_cells)):
            if nuclei_profile.has_specific_gene_cell_expression_with_feature(intron_genes[i],intron_cells[j],'intron'):
                combinations+=1

    matrix_outf = open(current_wd+'/intron_expression/matrix.mtx','w')
    matrix_outf.write('%%MatrixMarket matrix coordinate integer general\n%\n')
    matrix_outf.write(str(len(intron_genes))+' '+str(len(intron_cells))+' '+str(combinations)+'\n')
    for i in range(0,len(intron_genes)):
        for j in range(0,len(intron_cells)):
            if nuclei_profile.has_specific_gene_cell_expression_with_feature(intron_genes[i],intron_cells[j],'intron'):
                matrix_outf.write(str(i+1)+' '+str(j+1)+' '+str(len(nuclei_profile.get_gene_cell_expression_with_feature(intron_genes[i],intron_cells[j],'intron')))+'\n')
    matrix_outf.close()


if args.debug is not None:
    [exon_percentile, intron_percentile] = define_graphic_region("ENSMUSG00000051951")
    print('Exon : ', exon_percentile)
    print('Norm Exon : ', _normalize_coverage(exon_percentile))
    print('----------------------------------------------------------')
    print('Intron : ', intron_percentile)
    print('Norm Intron : ', _normalize_coverage(intron_percentile))

if args.coverageAnalysis is not None:
    top_genes = []
    exon_intron_df = pd.DataFrame()
    top_exonic_counts = sorted_gene_expression_by_exon[0:int(args.coverageAnalysis)]
    for i,v in top_exonic_counts:
        top_genes.append(i)
    print('Selcted top genes for Coverage analysis : ', len(top_genes))

    with multiprocessing.Pool(processes=args.core) as ca_p:
        ca_outputs = list(tqdm.tqdm(ca_p.imap(define_graphic_region,top_genes),total=len(top_genes)))
    ca_p.close()
    ca_p.join()
    ca_outputs = list(filter(lambda x: x is not False, ca_outputs))         #Remove exon/intron length < 300
    print('Remaining top genes for Coverage analysis : ', len(ca_outputs))
    all_exons_cov_percentile_combined = _combine_normalized_coverage_percentile(list(zip(*ca_outputs))[0])
    all_introns_cov_percentile_combined = _combine_normalized_coverage_percentile(list(zip(*ca_outputs))[1])

    norm_all_exons_cov_percentile_combined = _normalize_coverage(all_exons_cov_percentile_combined)
    norm_all_introns_cov_percentile_combined = _normalize_coverage(all_introns_cov_percentile_combined)

    #exon_intron_df = pd.DataFrame()


    ca_outf = open('top_selected_merged_genes_exon_intron_coverageAnalysis.txt','w')
    ca_outf.write('exons_coverage\t')
    for i in range(0,len(all_exons_cov_percentile_combined)):
        #exon_intron_df = exon_intron_df.append(pd.DataFrame.from_dict({i:['exon',all_exons_cov_percentile_combined[i]]}))
        ca_outf.write(str(all_exons_cov_percentile_combined[i])+'\t')
    ca_outf.write('\n')
    ca_outf.write('introns_coverage\t')
    for j in range(0,len(all_introns_cov_percentile_combined)):
        #exon_intron_df = exon_intron_df.append(pd.DataFrame.from_dict({i:['intron',all_introns_cov_percentile_combined[j]]}))
        ca_outf.write(str(all_introns_cov_percentile_combined[j])+'\t')
    ca_outf.write('\n')
    ca_outf.write('SCALED_exons_coverage\t')
    for i in range(0,len(norm_all_exons_cov_percentile_combined)):
        #exon_intron_df = exon_intron_df.append(pd.DataFrame.from_dict({i:['exon',all_exons_cov_percentile_combined[i]]}))
        ca_outf.write(str(norm_all_exons_cov_percentile_combined[i])+'\t')
    ca_outf.write('\n')
    ca_outf.write('SCALED_introns_coverage\t')
    for j in range(0,len(norm_all_introns_cov_percentile_combined)):
        #exon_intron_df = exon_intron_df.append(pd.DataFrame.from_dict({i:['intron',all_introns_cov_percentile_combined[j]]}))
        ca_outf.write(str(norm_all_introns_cov_percentile_combined[j])+'\t')
    ca_outf.write('\n')
    ca_outf.close()

    #exon_intron_df.columns = ['Percentile','Type','Normalized_coverage']

    x1 = pd.Series(all_exons_cov_percentile_combined,name='Exon coverage')
    x2 = pd.Series(all_introns_cov_percentile_combined,name='Intron coverage')

    if args.no is None:
        sns.set(style='white',palette='muted')
        fig, ax = plt.subplots()
        ax = sns.jointplot(x1,x2,kind='kde',size=7,space=0, ylim=(1,100),xlim=(1,100))
        #sns.tsplot(data=exon_intron_df, time="Percentile", value="Normalized_coverage",condition="type")
        fig.savefig('jplot.pdf')
        #sns.plt.show()
    #print('Exon data: ')
    #print(all_exons_cov_percentile_combined)
    #print('------------------------------------------------------')
    #print('Intron data: ')
    #print(all_introns_cov_percentile_combined)

    print('End of Coverage Analysis')

if args.express is not None:
    top_genes = []
    exon_intron_df = pd.DataFrame()
    top_exonic_counts = sorted_gene_expression_by_exon[0:int(args.express)]
    for i,v in top_exonic_counts:
        top_genes.append(i)
    print('Selcted top genes for intronic modeling : ', len(top_genes))

    with multiprocessing.Pool(processes=args.core) as ex_p:
        ex_outputs = list(tqdm.tqdm(ex_p.imap(intron_worker,top_genes),total=len(top_genes)))
    ex_p.close()
    ex_p.join()

    print('End of intronic profiling')

    ind = merging_exon_intron_expression_dict_CM(list(list(zip(*ex_outputs))[1]))
    top_recipricol_genes = sorted_gene_expression_by_exon[0:len(top_genes)]
    outf = open('intron_based_individual_expression.txt','w')
    outf.write('#Zero-based_coordinate_system\n')
    outf.write('gene_id\tgene_name\tstrand\tchromosome\texon_start\texon_end\texon_unique_read_count\n')
    for i,v in top_recipricol_genes:
        #outf.write(i+'\t'+mm10_gtf.get_gene_name(i)+'\t'+mm10_gtf.get_gene_strand(i)+'\t')
        intron_blocks = ind[i]
        for block in intron_blocks:
            outf.write(i+'\t'+mm10_gtf.get_gene_name(i)+'\t'+mm10_gtf.get_gene_strand(i)+'\t'+str(block[0])+'\t'+str(block[1])+'\t'+str(block[2])+'\t'+str(block[3])+'\n')
    outf.close()

    ##After sorting, I dont need to preserve order. Convert list of tuples to dictionary for the ease of writing
    sorted_gene_expression_by_exon = dict(sorted_gene_expression_by_exon)

    merged_top_gene_intron_expression = merging_exon_intron_expression_dict_CM(list(list(zip(*ex_outputs))[0]))
    #print(ex_outputs)
    outf = open('top_selected_genes_exon_intron_level.txt','w')
    outf.write('gene_id\tgene_name\texon_counts\tintron_counts\n')
    for tgene in top_genes:
        texon = sorted_gene_expression_by_exon[tgene]
        tintron = merged_top_gene_intron_expression[tgene]
        outf.write(tgene+'\t'+mm10_gtf.get_gene_name(tgene)+'\t'+str(texon)+'\t'+str(tintron)+'\n')
    outf.close()

    if args.no is None:
        for sel_gene in top_genes:
            exon_intron_df = exon_intron_df.append(pd.DataFrame.from_dict({sel_gene:[sorted_gene_expression_by_exon[sel_gene],merged_top_gene_intron_expression[sel_gene]]},orient='index'))
        #for curated_gene in list(list(zip(*ex_outputs))[1]):
        #    exon_intron_df = exon_intron_df.append(pd.DataFrame.from_dict(curated_gene,orient='index'))
        columns=['exon','intron']
        exon_intron_df.columns = columns
        fig, ax = plt.subplots()
        sns.set_style("whitegrid")
        ax = sns.jointplot(x='intron', y='exon', data=exon_intron_df,color='r',size=7)
        fig.savefig('jointplot.png')
        sns.plt.show()
