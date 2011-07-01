#!/usr/bin/env python
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2011 Vipin T. Sreedharan
# Copyright (C) 2009-2011 Friedrich Miescher Laboratory of the Max Planck Society
#
# Description: This program is used to get genome annotation from a valid GFF3 formated file.
# 

import re, sys, time
import numpy as np
import scipy.io

def addCDSphase(strand, cds):
    """Add CDS phase to the CDS exons"""
    
    cds_region, cds_flag = [], 0 
    if strand == '+':
        for cdspos in cds:
            if cds_flag == 0:
                cdspos = (cdspos[0], cdspos[1], 0)
                diff = (cdspos[1]-(cdspos[0]-1))%3
            else:
                xy = 0
                if diff == 0: 
                    cdspos = (cdspos[0], cdspos[1], 0)
                elif diff == 1: 
                    cdspos = (cdspos[0], cdspos[1], 2)
                    xy = 2
                elif diff == 2: 
                    cdspos = (cdspos[0], cdspos[1], 1)
                    xy = 1
                diff = ((cdspos[1]-(cdspos[0]-1))-xy)%3
            cds_region.append(cdspos)
            cds_flag = 1 
    elif strand == '-':
        cds.reverse()
        for cdspos in cds: 
            if cds_flag == 0:
                cdspos = (cdspos[0], cdspos[1], 0)
                diff = (cdspos[1]-(cdspos[0]-1))%3
            else:  
                xy = 0 
                if diff == 0: 
                    cdspos = (cdspos[0], cdspos[1], 0)
                elif diff == 1:
                    cdspos = (cdspos[0], cdspos[1], 2)
                    xy = 2
                elif diff == 2: 
                    cdspos = (cdspos[0], cdspos[1], 1)
                    xy = 1
                diff = ((cdspos[1]-(cdspos[0]-1))-xy)%3
            cds_region.append(cdspos)
            cds_flag = 1
        cds_region.reverse()
    return cds_region

def addUTR3Exon(utr_3, exon):
    """Add UTR3 region present in outside of exon coordinates and make complete exon at 3' end."""
    
    if exon[-1][1] == utr_3[0][0]:
        jun_cod = [exon[-1][0], utr_3[0][1]]
        exon[-1] = jun_cod
        utr_3 = utr_3[1:]
        for p in utr_3:exon.append(p)
    else:
        for p in utr_3:exon.append(p)
    return exon

def addUTR5Exon(utr_5, exon):
    """Add UTR5 region present in out side of exon coordinates and make complete at 5' end."""
    
    if utr_5[-1][1] == exon[0][0]:
        jun_cod = [utr_5[-1][0], exon[0][1]]
        exon[0] = jun_cod
        utr_5 = utr_5[:-1]
        for p in utr_5:exon.insert(0, p)
    else:
        for p in utr_5:exon.insert(0, p)
    return exon

def createEXONcoordinate(strand_p, five_p_utr, cds_cod, three_p_utr):
    """Create exon cordinates from UTR's and CDS region"""
     
    exon_pos = []
    if strand_p == '+':        
        utr5_start, utr5_end = 0, 0
        if five_p_utr != []:utr5_start = five_p_utr[-1][0];utr5_end = five_p_utr[-1][1]
        cds_5start, cds_5end = cds_cod[0][0], cds_cod[0][1]
        jun_exon = []
        if cds_5start-utr5_end == 0 or cds_5start-utr5_end == 1:jun_exon = [utr5_start, cds_5end]    
        if len(cds_cod) == 1:
            five_prime_flag = 0
            if jun_exon != []:
                five_p_utr = five_p_utr[:-1]
                five_prime_flag = 1
            for utr5 in five_p_utr:exon_pos.append(utr5)
            jun_exon = []
            utr3_start, utr3_end = 0, 0
            if three_p_utr != []: 
                utr3_start, utr3_end = three_p_utr[0][0], three_p_utr[0][1]
            if utr3_start-cds_5end == 0 or utr3_start-cds_5end == 1:jun_exon = [cds_5start, utr3_end]
            three_prime_flag = 0
            if jun_exon != []: 
                cds_cod = cds_cod[:-1]
                three_p_utr = three_p_utr[1:]
                three_prime_flag = 1
            if five_prime_flag == 1 and three_prime_flag == 1:exon_pos.append([utr5_start, utr3_end])
            if five_prime_flag == 1 and three_prime_flag == 0:
                exon_pos.append([utr5_start, cds_5end])
                cds_cod = cds_cod[:-1]
            if five_prime_flag == 0 and three_prime_flag == 1:exon_pos.append([cds_5start, utr3_end])
            for cds in cds_cod:exon_pos.append(cds)
            for utr3 in three_p_utr:exon_pos.append(utr3)
        else:    
            if jun_exon != []:
                five_p_utr = five_p_utr[:-1]
                cds_cod = cds_cod[1:]
            for utr5 in five_p_utr:exon_pos.append(utr5)
            exon_pos.append(jun_exon) if jun_exon != [] else ''
            jun_exon, utr3_start, utr3_end = [], 0, 0
            if three_p_utr != []:
                utr3_start = three_p_utr[0][0]
                utr3_end = three_p_utr[0][1]
            cds_3start, cds_3end = cds_cod[-1][0], cds_cod[-1][1]
            if utr3_start-cds_3end == 0 or utr3_start-cds_3end == 1:jun_exon = [cds_3start, utr3_end]
            if jun_exon != []:
                cds_cod = cds_cod[:-1]
                three_p_utr = three_p_utr[1:]
            for cds in cds_cod:exon_pos.append(cds)
            exon_pos.append(jun_exon) if jun_exon != [] else ''
            for utr3 in three_p_utr:exon_pos.append(utr3)
    elif strand_p == '-':
        utr3_start, utr3_end = 0, 0        
        if three_p_utr != []:
            utr3_start = three_p_utr[-1][0]
            utr3_end = three_p_utr[-1][1]
        cds_3start, cds_3end = cds_cod[0][0], cds_cod[0][1]
        jun_exon = []
        if cds_3start-utr3_end == 0 or cds_3start-utr3_end == 1:jun_exon = [utr3_start, cds_3end]  
        if len(cds_cod) == 1:    
            three_prime_flag = 0
            if jun_exon != []:
                three_p_utr = three_p_utr[:-1]
                three_prime_flag = 1
            for utr3 in three_p_utr:exon_pos.append(utr3)
            jun_exon, utr5_start, utr5_end = [], 0, 0
            if five_p_utr != []:
                utr5_start = five_p_utr[0][0]
                utr5_end = five_p_utr[0][1]
            if utr5_start-cds_3end == 0 or utr5_start-cds_3end == 1:jun_exon = [cds_3start, utr5_end]
            five_prime_flag = 0
            if jun_exon != []:
                cds_cod = cds_cod[:-1]
                five_p_utr = five_p_utr[1:]
                five_prime_flag = 1
            if three_prime_flag == 1 and five_prime_flag == 1:exon_pos.append([utr3_start, utr5_end])
            if three_prime_flag == 1 and five_prime_flag == 0:
                exon_pos.append([utr3_start, cds_3end])
                cds_cod = cds_cod[:-1]
            if three_prime_flag == 0 and five_prime_flag == 1:exon_pos.append([cds_3start, utr5_end])        
            for cds in cds_cod:exon_pos.append(cds)
            for utr5 in five_p_utr:exon_pos.append(utr5)
        else:
            if jun_exon != []:three_p_utr = three_p_utr[:-1];cds_cod = cds_cod[1:]
            for utr3 in three_p_utr:exon_pos.append(utr3)   
            if jun_exon != []:exon_pos.append(jun_exon)
            jun_exon, utr5_start, utr5_end = [], 0, 0
            if five_p_utr != []:
                utr5_start = five_p_utr[0][0]
                utr5_end = five_p_utr[0][1]    
            cds_5start, cds_5end = cds_cod[-1][0], cds_cod[-1][1]
            if utr5_start-cds_5end == 0 or utr5_start-cds_5end == 1:jun_exon = [cds_5start, utr5_end]
            if jun_exon != []:
                cds_cod = cds_cod[:-1]
                five_p_utr = five_p_utr[1:]
            for cds in cds_cod:exon_pos.append(cds)
            if jun_exon != []:exon_pos.append(jun_exon)    
            for utr5 in five_p_utr:exon_pos.append(utr5)
    return exon_pos

def init_gene():
    """Defining the gene information"""
    
    gene_info = dict(
          id = '',
          anno_id = [],
          confgenes_id = [],
          name = '',
          source = '',
          gene_info = {},
          alias = '',
          name2 = [],
          strand = '',
          chr = '',
          chr_num = [],
          paralogs = [],
          start = '',
          stop = '',
          transcripts = [],
          transcript_info = [],
          transcript_status = [],
          transcript_valid = [],
          exons = [],
          exons_confirmed = [],
          cds_exons = [],
          utr5_exons = [],
          utr3_exons = [],
          tis = [],
          tis_conf = [],
          tis_info = [],
          cdsStop = [],
          cdsStop_conf = [],
          cdsStop_info = [],
          tss = [],
          tss_info = [],
          tss_conf = [],
          cleave = [],
          cleave_info = [],
          cleave_conf = [],
          polya = [],
          polya_info = [],
          polya_conf = [],
          is_alt = [],
          is_alt_spliced = [],
          is_valid = [],
          transcript_complete = [],
          is_complete = [],
          is_correctly_gff3_referenced = '',
          splicegraph = []
          )
    return gene_info

def GFFParse(ref_file):
    """Extracting annotated features from a GFF file based on feature identifier mapping."""

    genes, transcripts, exons, utr5, utr3, cds = dict(), dict(), dict(), dict(), dict(), dict()
    ref_fh = open(ref_file, 'rU')
    for gln in ref_fh:
        gln = gln.strip('\n\r').split('\t')
        if not gln:continue ## not considering an empty line 
        if re.match(r'#', gln[0]) or re.match(r'>', gln[0]):continue ## not considering commented and FASTA header lines from GFF
        if len(gln) == 1:continue ## not considering if FASTA sequence along with GFF
        assert len(gln) == 9, '\t'.join(gln) ## a valid GFF line contains only 9 tab-delimited fields
        if '' in gln:continue ## empty fields in any line ? 
        # TODO 1: include all possible first level features
        if gln[2] == 'gene':
            gid, desc = None, dict(chr = gln[0], start = gln[3], stop = gln[4], orient = gln[6], src = gln[1])
            for atb in gln[-1].split(';'):
                if atb == '':continue
                atb = atb.split('=')
                if atb[0] == 'ID':gid = atb[1];continue
                desc[atb[0]] = atb[1]
            genes[(gln[0], gid)] = desc
        # TODO 2: include all possible second level features 
        elif gln[2] == 'mRNA':
            gid, desc = None, dict(chr = gln[0], start = gln[3], stop = gln[4], orient = gln[6], src = gln[1])
            for atb in gln[-1].split(';'):
                if atb == '':continue
                atb = atb.split('=')
                if atb[0] == 'Parent':gid = atb[1];continue
                desc[atb[0]] = atb[1]
            if (gln[0], gid) in transcripts:
                transcripts[(gln[0], gid)].append(desc)
            else:
                transcripts[(gln[0], gid)] = [desc]
        # TODO 3: get third level features
        elif gln[2] == 'exon':
            tid, desc = None, dict(chr = gln[0], start = gln[3], stop = gln[4], orient = gln[6], src = gln[1])
            for atb in gln[-1].split(';'):
                if atb == '':continue
                atb = atb.split('=')
                if atb[0] == 'Parent':tid = atb[1];continue
                desc[atb[0]] = atb[1]
            for fid in tid.split(','):
                if (gln[0], fid) in exons:
                    exons[(gln[0], fid)].append(desc)
                else:
                    exons[(gln[0], fid)] = [desc] 
        elif gln[2] == 'five_prime_UTR':
            tid, desc = None, dict(chr = gln[0], start = gln[3], stop = gln[4], orient = gln[6], src = gln[1])
            for atb in gln[-1].split(';'):
                if atb == '':continue
                atb = atb.split('=')
                if atb[0] == 'Parent':tid = atb[1];continue
                desc[atb[0]] = atb[1]
            for fid in tid.split(','):
                if (gln[0], fid) in utr5:
                    utr5[(gln[0], fid)].append(desc)
                else:
                    utr5[(gln[0], fid)] = [desc] 
        elif gln[2] == 'CDS':
            tid, desc = None, dict(chr = gln[0], start = gln[3], stop = gln[4], orient = gln[6], src = gln[1])
            for atb in gln[-1].split(';'):
                if atb == '':continue
                atb = atb.split('=')
                if atb[0] == 'Parent':tid = atb[1];continue
                desc[atb[0]] = atb[1]
            for fid in tid.split(','):
                if (gln[0], fid) in cds:
                    cds[(gln[0], fid)].append(desc)
                else:
                    cds[(gln[0], fid)] = [desc] 
        elif gln[2] == 'three_prime_UTR':
            tid, desc = None, dict(chr = gln[0], start = gln[3], stop = gln[4], orient = gln[6], src = gln[1])
            for atb in gln[-1].split(';'):
                if atb == '':continue
                atb = atb.split('=')
                if atb[0] == 'Parent':tid = atb[1];continue
                desc[atb[0]] = atb[1]
            for fid in tid.split(','):
                if (gln[0], fid) in utr3:
                    utr3[(gln[0], fid)].append(desc)
                else:
                    utr3[(gln[0], fid)] = [desc] 
    ref_fh.close()
    return genes, transcripts, exons, utr5, cds, utr3
    
def organizeGenes(genes, transcripts, exons, utr5, cds, utr3):
    """Generate an ordered gene list from GFF information"""
    
    gene_cnt, gene_models = 1, []
    for fid, fdet in sorted(genes.items()):
        if fid in transcripts:
            gene = init_gene()
            gene['id'] = gene_cnt
            gene['name'] = fid[1]
            gene['chr'] = fdet['chr']
            gene['source'] = fdet['src']
            if fdet['orient'] == '+': ## ?? check with GFF2Anno modules and GFF documentation
                gene['start'] = int(fdet['start'])
                gene['stop'] = int(fdet['stop']) + 1
            elif fdet['orient'] == '-':
                gene['start'] = int(fdet['start']) -1
                gene['stop'] = int(fdet['stop'])
            gene['strand'] = fdet['orient']
            if gene['strand'] != '+' and gene['strand'] != '-': gene['strand'] = '.'
            gn_info = dict()
            for key, val in fdet.items():
                if key in ['src', 'start', 'stop', 'chr', 'orient']:continue
                gn_info[key] = val
            if gn_info == {}:gn_info['ID'] = fid[1]
            gene['gene_info'] = gn_info
            if len(transcripts[fid]) > 1:
                gene['is_alt_spliced'] = 1
                gene['is_alt'] = 1
            else:
                gene['is_alt_spliced'] = 0
                gene['is_alt'] = 0
            for ftid in sorted(transcripts[fid]):
                gene['transcripts'].append(ftid['ID'])
                
                exon_cod, utr5_cod, utr3_cod, cds_reg = [], [], [], []

                if (gene['chr'], ftid['ID']) in exons:
                    if gene['strand'] == '+': ## ??
                        for ele_det in exons[(gene['chr'], ftid['ID'])]:exon_cod.append([int(ele_det['start']), int(ele_det['stop'])+1])
                    elif gene['strand'] == '-':
                        for ele_det in exons[(gene['chr'], ftid['ID'])]:exon_cod.append([int(ele_det['start'])-1, int(ele_det['stop'])])
                if (gene['chr'], ftid['ID']) in utr5: 
                    if gene['strand'] == '+':
                        for ele_det in utr5[(gene['chr'], ftid['ID'])]:utr5_cod.append([int(ele_det['start']), int(ele_det['stop'])+1])
                    elif gene['strand'] == '-':
                        for ele_det in utr5[(gene['chr'], ftid['ID'])]:utr5_cod.append([int(ele_det['start'])-1, int(ele_det['stop'])])
                if (gene['chr'], ftid['ID']) in utr3: 
                    if gene['strand'] == '+':
                        for ele_det in utr3[(gene['chr'], ftid['ID'])]:utr3_cod.append([int(ele_det['start']), int(ele_det['stop'])+1])
                    elif gene['strand'] == '-':
                        for ele_det in utr3[(gene['chr'], ftid['ID'])]:utr3_cod.append([int(ele_det['start'])-1, int(ele_det['stop'])])
                if (gene['chr'], ftid['ID']) in cds: 
                    if gene['strand'] == '+':
                        for ele_det in cds[(gene['chr'], ftid['ID'])]:cds_reg.append([int(ele_det['start']), int(ele_det['stop'])+1])
                    elif gene['strand'] == '-': ## ?? check with GFF2Anno modules.
                        for ele_det in cds[(gene['chr'], ftid['ID'])]:cds_reg.append([int(ele_det['start'])-1, int(ele_det['stop'])])
                ## make feature coordinates in generalized manner, GFF file may contain ascending or descending order, Here we are making in ASCENDING order.
                if gene['strand'] == '-':
                    if exon_cod != [] and len(exon_cod) != 1:
                        if exon_cod[0][0] > exon_cod[-1][0]: exon_cod.reverse()
                    if cds_reg != [] and len(cds_reg) != 1:
                        if cds_reg[0][0] > cds_reg[-1][0]: cds_reg.reverse()
                    if utr3_cod != [] and len(utr3_cod) != 1:
                        if utr3_cod[0][0] > utr3_cod[-1][0]: utr3_cod.reverse()
                    if utr5_cod != [] and len(utr5_cod) != 1:
                        if utr5_cod[0][0] > utr5_cod[-1][0]: utr5_cod.reverse()
                ## create exon coordiantes from UTR's and CDS region
                if exon_cod == []:
                    if cds_reg != []:exon_cod = createEXONcoordinate(gene['strand'], utr5_cod, cds_reg, utr3_cod)
                ## create complete exon coordinates if UTR's are splitted TODO check for a code clean 
                if gene['strand'] == '+':
                    if utr5_cod != []:
                        if utr5_cod[-1][1] <= exon_cod[0][0]: exon_cod = addUTR5Exon(utr5_cod, exon_cod)
                    if utr3_cod != []:
                        if exon_cod[-1][1] <= utr3_cod[0][0]: exon_cod = addUTR3Exon(utr3_cod, exon_cod)
                elif gene['strand'] == '-':
                    if utr5_cod != []:
                        if exon_cod[-1][1] <= utr5_cod[0][0]: exon_cod = addUTR3Exon(utr5_cod, exon_cod)
                    if utr3_cod != []:
                        if utr3_cod[-1][1] <= exon_cod[0][0]: exon_cod = addUTR5Exon(utr3_cod, exon_cod)
                ## find [tis, cdsStop] and [tss, cleave] for the transcript
                tis, cdsStop, tss, cleave = [], [], [], []
                if cds_reg != []:
                    if gene['strand'] == '+':
                        tis = [cds_reg[0][0]]
                        cdsStop = [cds_reg[-1][1]-3]
                    elif gene['strand'] == '-':
                        tis = [cds_reg[-1][1]]
                        cdsStop = [cds_reg[0][0]+3]
                if utr5_cod != []:
                    if gene['strand'] == '+':
                        tss = [utr5_cod[0][0]]
                    elif gene['strand'] == '-':
                        tss = [utr5_cod[-1][1]]
                if utr3_cod != []: 
                    if gene['strand'] == '+':
                        cleave = [utr3_cod[-1][1]]
                    elif gene['strand'] == '-':
                        cleave = [utr3_cod[0][0]]
                ## calculating CDS phase 
                if cds_reg != []:cds_cod_phase = addCDSphase(gene['strand'], cds_reg)
                ## transctipt status 
                cds_status, exon_status, utr_status = 0, 0, 0
                if cds_cod_phase != []: cds_status = 1
                if exon_cod != []: exon_status = 1
                if utr5_cod != [] or utr3_cod != []: utr_status = 1
                if cds_status != 0 and exon_status != 0 and utr_status != 0: 
                    gene['transcript_status'].append(1)
                else: 
                    gene['transcript_status'].append(0)
                gene["exons"].append(exon_cod)
                gene['cds_exons'].append(cds_cod_phase)
                gene['utr3_exons'].append(utr3_cod)
                gene['utr5_exons'].append(utr5_cod)
                gene['tis'].append(tis)
                gene['cdsStop'].append(cdsStop)
                gene['tss'].append(tss)
                gene['cleave'].append(cleave)
                #break ## one transcript 
            gene['is_correctly_gff3_referenced'] = 1 # TODO gff3 id mapping
            gene = processFeature(gene)
            gene_models.append(gene)
            gene_cnt += 1
        #break ## one gene 
    return gene_models
        
def processFeature(singlegene):
    """Make feature value compactable to write in a .mat format"""

    comp_exon = np.zeros((len(singlegene['exons']),), dtype=np.object)
    for i in range(len(singlegene['exons'])):
        comp_exon[i]= np.array(singlegene['exons'][i])
    singlegene['exons'] = comp_exon
    comp_cds = np.zeros((len(singlegene['cds_exons']),), dtype=np.object)
    for i in range(len(singlegene['cds_exons'])):
        comp_cds[i]= np.array(singlegene['cds_exons'][i])
    singlegene['cds_exons'] = comp_cds
    comp_utr3 = np.zeros((len(singlegene['utr3_exons']),), dtype=np.object)
    for i in range(len(singlegene['utr3_exons'])):
        comp_utr3[i]= np.array(singlegene['utr3_exons'][i])
    singlegene['utr3_exons'] = comp_utr3
    comp_utr5 = np.zeros((len(singlegene['utr5_exons']),), dtype=np.object)
    for i in range(len(singlegene['utr5_exons'])):
        comp_utr5[i]= np.array(singlegene['utr5_exons'][i])
    singlegene['utr5_exons'] = comp_utr5
    comp_transcripts = np.zeros((len(singlegene['transcripts']),), dtype=np.object)
    for i in range(len(singlegene['transcripts'])):
        comp_transcripts[i]= np.array(singlegene['transcripts'][i])
    singlegene['transcripts'] = comp_transcripts
    comp_tss = np.zeros((len(singlegene['tss']),), dtype=np.object)
    for i in range(len(singlegene['tss'])):
        comp_tss[i]= np.array(singlegene['tss'][i])
    singlegene['tss'] = comp_tss
    comp_tis = np.zeros((len(singlegene['tis']),), dtype=np.object)
    for i in range(len(singlegene['tis'])):
        comp_tis[i]= np.array(singlegene['tis'][i])
    singlegene['tis'] = comp_tis
    comp_cleave = np.zeros((len(singlegene['cleave']),), dtype=np.object)
    for i in range(len(singlegene['cleave'])):
        comp_cleave[i]= np.array(singlegene['cleave'][i])
    singlegene['cleave'] = comp_cleave
    comp_cdsStop = np.zeros((len(singlegene['cdsStop']),), dtype=np.object)
    for i in range(len(singlegene['cdsStop'])):
        comp_cdsStop[i]= np.array(singlegene['cdsStop'][i])
    singlegene['cdsStop'] = comp_cdsStop
    return singlegene

def __main__():
    
    try:
        gff_file = sys.argv[1]
        mat_file = sys.argv[2]
    except:
        sys.stderr.write('ERROR: Provide a valid GFF/GTF file and result file name in .mat extension\nUSAGE: GFFParser.py <GFF file> <gene_list.mat>\n')
        sys.exit(-1)
    #print time.asctime( time.localtime(time.time()) )
    ## extract genome annotations from GFF file 
    genes, transcripts, exons, utr5, cds, utr3 = GFFParse(gff_file)
    ## generate a Gene model structure 
    genes_data = organizeGenes(genes, transcripts, exons, utr5, cds, utr3)
    ## saving gene models into .mat file 
    scipy.io.savemat(mat_file, mdict={'gene_models':genes_data}, format='5', oned_as='row')
    #print time.asctime( time.localtime(time.time()) )

if __name__=="__main__":__main__()
