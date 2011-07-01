## File name : GFFFastParser.py
## Description : Using Python GFF parser, parse the huge GFF3 file 

## Modules used
import GFFParser

def create_dict_from_feature(feature,gene_cont):

    gene = init_gene()
    gene["id"] = gene_count
    gene["name"] = feature.id    
    # gene[alias] = 
    if feature.strand == 1:
       sign = '+'
    elif feature.strand == -1:
       sign = '-'  
    gene["strand"] = sign
    gene["chr"] = final.id
    gene["chr_num"] = final.id
    gene["start"] = feature.location._start.position
    gene["stop"] = feature.location._end.position
    
    for sub in feature.sub_features:
        gene["transcripts"] = sub.id 
        # gene["transcript_status"] = 
        # gene["transcript_valid"] =
        
        rows = len(sub.sub_features)
        exon_pos = [[0 for i in range(2)] for j in range(rows)]
        i = 0
        for child in sub.sub_features:
            e_start = child.location._start.position
            e_end = child.location._end.position
            for e_s in range(rows):
                for e_e in range(2):
                    if i == e_s:
                        if e_e == 0:
                            exon_pos[e_s][e_e] = e_start
                    if i == e_s:
                        if e_e == 1:            
                            exon_pos[e_s][e_e] = e_end
            i += 1
        
        gene["exons"] = exon_pos  
        # gene["cds_exons"] =
        # gene["utr5_exons"] =  
        # gene["utr3_exons"] =        
        # gene["tis"] =  
        # gene["cdsStop"] =
        # gene["tss"] =  
        # gene["cleave"] =   
    
    return gene

def init_gene():

    gfffields = dict(
          id = [],
          anno_id = [],
          confgenes_id = [],
          name = '',
          alias = [],
          name2 = [],
          strand = '',
          chr = '',
          chr_num = [],
          paralogs = [],
          start = [],
          stop = [],
          transcripts = {},
          transcript_info = [],
          transcript_status = [],
          transcript_valid = [],
          exons = {},
          exons_confirmed = {},
          cds_exons = {},
          utr5_exons = {},
          utr3_exons = {},
          tis = [],
          tis_conf = [],
          tis_info = [],
          cdsStop = [],
          cdsStop_conf = [],
          cdsStop_info = [],
          transacc = [],
          transacc_info = [],
          transacc_conf = [],
          tss = [],
          tss_info = [],
          tss_conf = [],
          cleave = [],
          cleave_info = [],
          cleave_conf =[],
          polya = [],
          polya_info = [],
          polya_conf = [],
          is_alt = [],
          is_alt_spliced = [],
          is_valid = [],
          is_complete = [],
          is_correctly_gff3_referenced = [],
          splicegraph = [],
	  in_operon = []
          )
    return gfffields
 
pgff = GFFParser.GFFMapReduceFeatureAdder(dict(), None)

exon_limit_info = dict(
	gff_type = ["gene","mRNA","exon"],
	gff_id = ["II"]
	)

pgff.add_features('../PythonGFF/small.gff3', exon_limit_info)

pgff.base["II"]

final = pgff.base["II"]

genes = []

gene_count = 1

for feature in final.features:
    
    gene = create_dict_from_feature(feature,gene_count)
    print gene["id"]
    print gene["name"]  
    print gene["strand"]  
    print gene["chr"]
    print gene["chr_num"]
    print gene["start"] 
    print gene["stop"]  
    print gene["transcripts"]
    print gene["exons"]  
    gene_count+=1
    break

#   genes.append(gene)

# save genes as mat file
