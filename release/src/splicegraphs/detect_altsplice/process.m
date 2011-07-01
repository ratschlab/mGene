addpath ../../utils
addpath ../../splicegraphs
addpath /fml/ag-raetsch/share/software/matlab_tools/utils

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('initializing genome\n') ;
genome_info = init_genome('/fml/ag-raetsch/share/projects/splicing/elegans_WS147_new/genome.config') ;
%genome_info = init_genome('/fml/ag-raetsch/share/projects/splicing/arabidopsis_new/genome.config') ;
%genome_info = init_genome('/fml/ag-raetsch/share/projects/splicing/drosophila_new/genome.config') ;
%genome_info = init_genome('/fml/ag-raetsch/share/projects/splicing/rerio_new/genome.config') ;
fprintf('directory: %s\n',genome_info.basedir) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CONFIRMED SEQUENCES %%%
%load(sprintf('%s/confirmed_sequences_merge.mat', genome_info.basedir)) ;
load(sprintf('%s/confirmed_sequences_infer.mat', genome_info.basedir)) ;


%%% CONFIRMED AND ANNOTATED SEQUENCES %%%
%load(sprintf('%s/conf_anno_sequences_merge.mat', genome_info.basedir)) ;
%load(sprintf('%s/conf_anno_sequences_infer.mat', genome_info.basedir)) ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% what has to be done %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1
  disp('Sorting exons...');
  for ix=1:length(genes)
    [dummy,exon_order] = sort(genes(ix).splicegraph{1}(1,:),2,'ascend');
    genes(ix).splicegraph{1} = genes(ix).splicegraph{1}(:,exon_order);
    genes(ix).splicegraph{2} = genes(ix).splicegraph{2}(exon_order,exon_order);
  end
end

%%% detect alternatively spliced genes %%%

%%% VERSION I  %%%
%[sum_alt, sum_const, idx_alt, idx_con, genes] = detect_altsplicing(genes) ;

%%% VERSION II %%%
genes = alt_const(genes) ;




% which rules hold %
if 0,
  detect_rule(genes, idx_alt) ;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% what can be be done %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% count new exons created by infer %%%
if 1
  [genes, new_exon_count]=detect_newexons(genes, idx_alt) ;
end


% detect exon skips %
if 1
  [idx_exon_skips, exon_exon_skips] = detect_exonskips(genes, idx_alt);
  %genes_exon_skips = genes(idx_exon_skips);
  %save(sprintf('%s/counting/genes.exon_skips.graph.mat',genome_info.basedir),'genes_exon_skips', 'exon_exon_skips') ;
end
  

%%% detect alternative 5 and 3 prime sites %%%
if 1
  [idx_alt_5prime,exon_alt_5prime, idx_alt_3prime,exon_alt_3prime] = detect_altprime(genes, idx_alt);
  genes_alt_5prime= genes(idx_alt_5prime);
  %save(sprintf('%s/counting/genes.alt_5prime.graph.mat',genome_info.basedir), 'genes_alt_5prime', 'exon_alt_5prime') ;
  genes_alt_3prime = genes(idx_alt_3prime);
  %save(sprintf('%s/counting/genes.alt_3prime.graph.mat',genome_info.basedir), 'genes_alt_3prime', 'exon_alt_3prime') ;
end


%%% detect intron retentions %%%
if 1
  [idx_intron_reten,intron_intron_reten] = detect_intronreten(genes, idx_alt) ;
  %genes_intron_reten = genes(idx_intron_reten);
  %save(sprintf('%s/counting/genes.intron_reten.graph.mat',genome_info.basedir), 'genes_intron_reten', 'exon_intron_reten') ;
end


%%% detect XOR exons %%%
if 1
  [idx_xor_exons, exon_xor_exons] = detect_xorexons(genes, idx_alt) ;
  %genes_xor_exons= genes(idx_xor_exons);
  %save(sprintf('%s/counting/genes.xor_exons.graph.mat',genome_info.basedir), 'genes_xor_exons', 'exon_xor_exons') ;
end


idx_multiple_skips=[] ;
%%% detect multiple exon skips %%%
if 1
  [idx_multiple_skips, exon_multiple_skips] = detect_multipleskips(genes, idx_alt) ;
  %genes_multiple_skips = genes(idx_multiple_skips);
  %save(sprintf('%s/counting/genes.multiple_skips.graph.mat',genome_info.basedir), 'genes_multiple_skips', 'exon_multiple_skips') ;
end


%%% detect incomplete exons in intronic regions %%%
if 1
  [idx_exon_intron5, exon_exon_intron5, idx_exon_intron3, exon_exon_intron3] = detect_exonintron(genes, idx_alt) ;
  %genes_exon_intron = genes(idx_exon_intron); 
  %save(sprintf('%s/counting/genes.exon_intron.graph.mat',genome_info.basedir), 'genes_exon_intron', 'exon_exon_intron') ;  
end


%%% detect alternative introns %%%
if 1
  [idx_alt_intron, introns_alt_intron] = detect_altintrons(genes, [1:length(genes)]) ;
  %genes_alt_intron = genes(idx_alt_intron); 
  %save(sprintf('%s/counting/genes.alt_intron.graph.mat',genome_info.basedir), 'genes_alt_intron', 'intron_alt_intron') ;
end


%%% detect undetermined splicing events %%%            
if 1
  idx_unknown = [];
  for ix=idx_alt
    if (mod(ix,50)==0)
      fprintf(1,'.');
    end
    if (isempty(find(ix==idx_exon_skips))) && ...
       (isempty(find(ix==idx_alt_5prime))) && ...
       (isempty(find(ix==idx_alt_3prime))) && ...
       (isempty(find(ix==idx_intron_reten))) && ...
       (isempty(find(ix==idx_alt_intron))) && ...  
       (isempty(find(ix==idx_xor_exons))) && ...
       (isempty(find(ix==idx_multiple_skips))) && ...
       (isempty(find(ix==idx_exon_intron5))) && ...
       (isempty(find(ix==idx_exon_intron3))),
      idx_unknown = [idx_unknown,ix];
    end
  end
  fprintf(1,'\n\nNumber of undetermined splicing events:\t\t\t\t%d\n', length(idx_unknown));
  %genes_unknown = genes(idx_unknown);
  %save(sprintf('%s/counting/genes.unknown.graph.mat',genome_info.basedir), 'genes_unknown') ;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% saving it the right way %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
  genes_detect={} ;
  for i=1:length(genes)
    if (mod(i,50)==0)
      fprintf(1,'.');
    end
    genes_detect(i).idx=i ;
    genes_detect(i).name=genes(i).name ;
    genes_detect(i).is_alt=genes(i).is_alt ;
    genes_detect(i).exonskips={} ;
    genes_detect(i).alt_5prime={} ;
    genes_detect(i).alt_3prime={} ;
    genes_detect(i).intron_reten={} ;
    genes_detect(i).alt_intron={} ;
    genes_detect(i).xor_exons={} ;
    genes_detect(i).incompl_intronreten5={} ;
    genes_detect(i).incompl_intronreten3={} ;
  end
  
  % exonskips % 
  k=1 ;
  for i=idx_exon_skips,
    genes_detect(i).exonskips(end+1).threeprimesite=genes(i).splicegraph{1}(2,exon_exon_skips(1,k)) ;
    genes_detect(i).exonskips(end).exon=[genes(i).splicegraph{1}(1, exon_exon_skips(2,k)), ...
                                         genes(i).splicegraph{1}(2, exon_exon_skips(2,k))] ;
    genes_detect(i).exonskips(end).fiveprimesite=genes(i).splicegraph{1}(1,exon_exon_skips(3,k)) ; 
    k=k+1 ;
  end
  
  % alternative 5 prime sites % 
  k=1 ;
  for i=idx_alt_5prime,
    if genes(i).strands(1)=='+',
      genes_detect(i).alt_5prime(end+1).threeprimesite=genes(i).splicegraph{1}(2,exon_alt_5prime(k).threeprimesite) ;
      genes_detect(i).alt_5prime(end).fiveprimesites=genes(i).splicegraph{1}(1, exon_alt_5prime(k).fiveprimesites) ;
    elseif genes(i).strands(1)=='-',
      genes_detect(i).alt_5prime(end+1).threeprimesite=genes(i).splicegraph{1}(1,exon_alt_5prime(k).threeprimesite) ;
      genes_detect(i).alt_5prime(end).fiveprimesites=genes(i).splicegraph{1}(2, exon_alt_5prime(k).fiveprimesites) ;
    end
    k=k+1 ;
  end
  
  % alternative 3 prime sites %
  k=1 ;
  for i=idx_alt_3prime,
    if genes(i).strands(1)=='+',
      genes_detect(i).alt_3prime(end+1).fiveprimesite=genes(i).splicegraph{1}(1,exon_alt_3prime(k).fiveprimesite) ;
      genes_detect(i).alt_3prime(end).threeprimesites=genes(i).splicegraph{1}(2,exon_alt_3prime(k).threeprimesites) ;
    elseif genes(i).strands(1)=='-',
      genes_detect(i).alt_3prime(end+1).fiveprimesite=genes(i).splicegraph{1}(2,exon_alt_3prime(k).fiveprimesite) ;
      genes_detect(i).alt_3prime(end).threeprimesites=genes(i).splicegraph{1}(1,exon_alt_3prime(k).threeprimesites) ;
    end
    k=k+1 ;
  end
  
  % detect intron retentions %
  k=1 ;
  for i=idx_intron_reten,
    genes_detect(i).intron_reten(end+1).intron=[genes(i).splicegraph{1}(2,intron_intron_reten(1,k)), ...
                                                genes(i).splicegraph{1}(1,intron_intron_reten(2,k))] ;
    k=k+1 ;
  end
  
  
  % detect alternative introns %
  k=1 ;
  for i=idx_alt_intron,
    genes_detect(i).alt_intron(end+1).intron1=[genes(i).splicegraph{1}(2,introns_alt_intron(1,k)), ...
                                               genes(i).splicegraph{1}(1,introns_alt_intron(2,k))] ;
    genes_detect(i).alt_intron(end).intron2=[genes(i).splicegraph{1}(2,introns_alt_intron(3,k)), ...
                                             genes(i).splicegraph{1}(1,introns_alt_intron(4,k))] ;    
    k=k+1 ;
  end
  
  % detect XOR exons %
  k=1 ;
  for i=idx_xor_exons,   
    genes_detect(i).xor_exons(end+1).threeprimesite=genes(i).splicegraph{1}(2,exon_xor_exons(1,k)) ;
    genes_detect(i).xor_exons(end).exon1=[genes(i).splicegraph{1}(1,exon_xor_exons(2,k)), ...
                                          genes(i).splicegraph{1}(2,exon_xor_exons(2,k))] ;
    genes_detect(i).xor_exons(end).exon2=[genes(i).splicegraph{1}(1,exon_xor_exons(3,k)), ...
                                          genes(i).splicegraph{1}(2,exon_xor_exons(3,k))] ;
    genes_detect(i).xor_exons(end).fiveprimesite=genes(i).splicegraph{1}(1,exon_xor_exons(4,k)) ;    
    k=k+1 ;
  end
  
  % detect incomplete exons in intronic regions (5 prime site) %
  k=1 ;
  for i=idx_exon_intron5,
    genes_detect(i).incompl_intronreten5(end+1).intronstart=genes(i).splicegraph{1}(2,exon_exon_intron5(1,k)) ;
    genes_detect(i).incompl_intronreten5(end).intronend=genes(i).splicegraph{1}(1,exon_exon_intron5(2,k)) ;
    genes_detect(i).incompl_intronreten5(end).exon=[genes(i).splicegraph{1}(1,exon_exon_intron5(3,k)), ...
                                                    genes(i).splicegraph{1}(2,exon_exon_intron5(3,k))] ;   
    k=k+1 ;
  end
  
  % detect incomplete exons in intronic regions (3 prime site) %
  k=1 ;
  for i=idx_exon_intron3,
    genes_detect(i).incompl_intronreten3(end+1).intronstart=genes(i).splicegraph{1}(2,exon_exon_intron3(1,k)) ;
    genes_detect(i).incompl_intronreten3(end).intronend=genes(i).splicegraph{1}(1,exon_exon_intron3(2,k)) ;
    genes_detect(i).incompl_intronreten3(end).exon=[genes(i).splicegraph{1}(1,exon_exon_intron3(3,k)), ...
                                                    genes(i).splicegraph{1}(2,exon_exon_intron3(3,k))] ;
    k=k+1 ;
  end
  
  save(sprintf('%s/genes_detect.mat',genome_info.basedir), 'genes_detect') ;
  
end



