clear;
fprintf(1,'Reading .gff files and constructing list of genes\n');
addpath ../sgetools/
read_gff_sge;
%read_gff_files;

%fprintf(1,'Loading Data');
%load Data/genes.all.mat

fprintf(1,'Remove ESTs with only one exon, and add chromosome and strand to name\n');
filter_single_exons;

fprintf(1,'Merge colocated genes\n');
gen_len = 0;
while gen_len ~= length(genes)
  gen_len = length(genes);
  merge_transcripts;  % merge colocated genes
end

save Data/elegans.genes.all.mat genes




%fprintf(1,'Include cDNA information\n');
%genes = check_genes(genes);
%save Data/elegans.genes.all.checked.mat genes







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We construct the splice graphs and merge the information. We
% detect some simple alternative splicing events, and take note on
% which genes that the occur. Also keep track of which
% corresponding exons that are involved in these events.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'Build splice graph\n');
build_splice_graph;
grow_splice_graph;
compact_splice_graph;
countspliceform;
%viewsplicegraph(genes(1));

save Data/elegans.genes.graph.mat genes idx_alt idx_alt_3prime ...
    idx_alt_5prime idx_con idx_exon_intron idx_exon_skips idx_intron_reten ...
    idx_multiple_skips idx_unknown idx_xor_exons ...
    exon_alt_3prime exon_alt_5prime exon_exon_intron exon_exon_skips ...
    exon_intron_reten exon_multiple_skips exon_xor_exons








%load genes.alt.est.mat
%allgenes = genes;

%clear genes;
%load genes.alt.mat
%allgenes = [allgenes,genes];

%clear genes;
%load genes.noalt.est.mat
%allgenes = [allgenes,genes];

%clear genes;
%load genes.noalt.mat
%allgenes = [allgenes,genes];

%clear genes;
%genes = allgenes;
%clear allgenes;


% account for duplicated gene descriptors

%[tmp,tmp,idx]=unique({genes(:).name});
%[list,idx2] = sort(idx) ;
%take_map = ones(1,length(genes)) ;
%for i=1:length(idx)-1,
%  %if list(i)==list(i+2), keyboard ; end ;
%  j=0 ; happened = 1 ;
%  while happened,
%    happened = 0 ;
%    j=j+1 ;
%    if list(i)==list(i+j), 
%      %keyboard
%      happened = 1 ;
%      assert(genes(idx2(i)).strands(1)==genes(idx2(i+j)).strands(1)) ;
%      assert(equal(genes(idx2(i)).chr, genes(idx2(i+j)).chr)) ;
%      genes(idx2(i)).is_alt = genes(idx2(i)).is_alt | genes(idx2(i+j)).is_alt;
%      genes(idx2(i)).transcripts = {genes(idx2(i)).transcripts{:} genes(idx2(i+j)).transcripts{:}} ;
%      genes(idx2(i)).exons = {genes(idx2(i)).exons{:} genes(idx2(i+j)).exons{:}} ;
%      genes(idx2(i)).strands = [genes(idx2(i)).strands genes(idx2(i+j)).strands] ;
%      take_map(idx2(i+j)) = 0 ;
%    end ;
%  end ;
%end ;
%idx = find(take_map) ;
%genes=genes(idx) ;

%save genes.all.mat

%build_splice_graph;
%grow_splice_graph;
%compact_splice_graph;

%save genes.all.graph.mat
