function genes = init_genes(num_genes,PAR,genes)
% genes = init_genes(num_genes,PAR,genes)
  
if nargin<1
  num_genes=0;
end
if nargin<2 || isempty(PAR),
  PAR.organism.name =[];
end
    
genes.id = [];
genes.anno_id = [];
genes.confgenes_id = [];
genes.name = '' ;
genes.alias = [];
genes.name2 = []; 
genes.strand = '';
genes.chr = '';
genes.chr_num = [];
genes.paralogs = [];

genes.start = [];
genes.stop = [];

genes.transcripts = {};
genes.transcript_info = [];
genes.transcript_status = [];
genes.transcript_valid = [];
genes.exons = {}; 
genes.exons_confirmed = {};
genes.cds_exons = {}; 
genes.utr5_exons = {}; 
genes.utr3_exons = {}; 
genes.tis = {};
genes.tis_conf = [];
genes.tis_info = [];

genes.cdsStop = {};
genes.cdsStop_conf = [];
genes.cdsStop_info = [];

if isequal(PAR.organism.name,'C_elegans')||(isfield(PAR.organism,'clade')&&isequal(PAR.organism.clade,'nematode'))
  genes.transacc = {};
  genes.transacc_info = [];
  genes.transacc_conf = [];
end
genes.tss = {};
genes.tss_info = [];
genes.tss_conf = [];
genes.cleave = {};
genes.cleave_info = [];
genes.cleave_conf = [];
genes.polya = {};
genes.polya_info = [];
genes.polya_conf = [];

genes.is_alt = [];
genes.is_alt_spliced = [];
genes.is_valid = [];
genes.transcript_complete = [];
genes.is_complete = [];
genes.is_correctly_gff3_referenced = [];

genes.splicegraph =[];

if isequal(PAR.organism.name,'C_elegans')||(isfield(PAR.organism,'clade')&&isequal(PAR.organism.clade,'nematode'))
  genes.in_operon=[];
end

if num_genes==0
  genes=genes([]);
elseif num_genes>1,
  genes(num_genes) = genes ;
end
