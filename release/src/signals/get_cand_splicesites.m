function [POS,LABEL,CONF] = get_cand_splicesites(P)
  
% [pos,label,Conf] = get_cand_splicesites(genes,seq,PAR)
% WARNING IMPLENENTATION OF COV_SKIPPED NOT FINISHED!!!

debug = 0;  
  
genes = P.genes ;
seq = P.seq ;
Signal = P.Signal ;
info = P.info ;
strand = P.strand ;
%clear P
  
% info = PAR.info_genes.(signal_name);
[info_names,sig,alt_names] = initialize_cand_search(Signal,info,genes,strand,'transcript',debug);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET ALL CONSENSUS SITES / INITIALIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

POS = find_consensus(seq,Signal,strand);

ALL = initialize_allconfs(POS, info_names);
ALL.LABEL =  -ones(1,length(POS)); 
clear POS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET TRANSSPLICE SPECIFIC STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isequal(Signal.name,'transacc')
  ALL = get_cand_transacc(genes,ALL,Signal,info_tran);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detecting exons in regions of alternative splicing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('detecting exons in regions of alternative splicing') ;
for i=1:length(genes)
  genes(i).alt_regions = detectaltregions(genes(i)) ;
  map = zeros(1,size(genes(i).splicegraph{1},2)) ;
  for j=1:size(genes(i).alt_regions,2),
    map = map | segment_overlap(genes(i).alt_regions(1,j)+1, ...
                                genes(i).alt_regions(2,j)-1, ...
                                genes(i).splicegraph{1}(1,:), ...
                                genes(i).splicegraph{1}(2,:)) ;
  end ;
  if ~isempty(genes(i).alt_regions),
    assert(any(map==1)) ;
  end ;
  
  %extent alt_region to neigbouring exons
  map = [0 map 0];
  idx_temp = find(map);
  map(idx_temp+1)=1;
  map(idx_temp-1)=1;
  genes(i).exon_in_alt_region=map(2:end-1) ;
  %% Version 4225 works without extension
end ;

genes = annotate_splicegraph(genes);
% genes = cover_splicegraph(genes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START LOOP OVER GENES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_invalid_genes = 0;
invalid_true = 0;

for g_idx = 1:length(genes)
  
  
  %fprintf(1,'%i\r',g_idx);
  gene = genes(g_idx) ;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % GET BASIC STUFF, 
  % if candidate pos is in gene, get gene id, if gene is valid, if
  % candidate pos is in alternative region 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [ALL,is_invalid] = basic_cand_stuff(gene,ALL,sig,debug);
  num_invalid_genes  = num_invalid_genes + is_invalid;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % GET TRUE SPLICE SITES
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [don_exons, acc_exons] = ind2sub(size(gene.splicegraph{2}), find(triu(gene.splicegraph{2}))) ;
  true_splice = [];
  
  if (isequal(sig,'acc') && strand=='+')||...
        (isequal(sig,'don') && strand=='-'),
    true_splice = unique(gene.splicegraph{1}(1,acc_exons)) ;
  elseif (isequal(sig,'don')&& strand=='+') || ...
        (isequal(sig,'acc')&& strand=='-')
    true_splice = unique(gene.splicegraph{1}(2,don_exons)) ;
  end
   
  [tmp,ii1,ii2] = intersect(ALL.POS,true_splice);
  
  if length(tmp)<length(true_splice)
    invalid_true = invalid_true+length(true_splice)-length(tmp) ;
  end
  if ~isequal(Signal.name,'transacc')
    ALL.LABEL(ii1) = 1;
  else
    ALL.ISSPLICE(ii1) = 1;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % GET TRUE ALTERNATIVE SPLICE SITES
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  if ~isempty(alt_names)
    [ALL,label_2] = get_true_alternative_sites(ALL,gene,true_splice,sig,debug);
  else
    label_2 = -ones(3,length(true_splice)); 
  end 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %%% CHECK COVERAGE 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  if ~isempty(info_names) 
    ALL = get_coverage(ALL, gene, true_splice, label_2, info, info_names);
  end
      
  % keyboard
end %%%loop over genes

if isequal(Signal.name,'transacc')
  ALL.COV1 = ALL.CONF_TRANS;
end
if nargout>=3,
  CONF = set_CONF(Signal, ALL, sig, info_names, '');
end ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLE DOWN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos=ALL.POS(ALL.LABEL==1);
neg_pos=ALL.POS(ALL.LABEL==-1);

% pos contains positive examples, neg_pos contains negative examples

%if Signal.sampling.downsample_neg~=1
%  fprintf(1,'Sampling down negatives:%d',Signal.sampling.downsample_neg);
%  randn('state',24593854);
%  dist = abs(ceil(randn(length(neg_pos),1)./Signal.sampling.downsample_neg/5+ Signal.sampling.mindist+1/Signal.sampling.downsample_neg));
%  dist(dist<Signal.sampling.mindist)=dist(dist<Signal.sampling.mindist)+Signal.sampling.mindist;
%  dd = cumsum(dist);
%  neg_pos = intersect(neg_pos,dd);
%elseif Signal.sampling.downsample_neg==1 & Signal.sampling.mindist>1
%  remove = 1;
%  while remove
%    dist = neg_pos(2:end)-neg_pos(1:end-1)
%    idx = find(dist<Signal.sampling.mindist)+1;
%    neg_pos(idx) = [];
%    if isempty(idx)
%      remove = 0;
%    end
%  end
%end

POS=[pos neg_pos];
LABEL=[ones(1,length(pos)), -ones(1,length(neg_pos))];

% keyboard
% assert(all(sum(Conf(1:3,LABEL==1),1)))
% assert(size(Conf,1)==length(Signal.Conf_names));
fprintf(1,'number of invalid genes: \t\t\t%i\n',num_invalid_genes)
fprintf(1,'number of true %s sites: \t\t\t%i\n',Signal.name,sum(LABEL==1))
fprintf(1,'number of invalid true sites (no consensus): \t%i\n',invalid_true)
fprintf(1,'number of decoy sites:\t\t\t\t%i\n',length(LABEL))
% fprintf(1,'mean confidence:\t\t\t\t%.2f\n',mean(sum(Conf,1)));

if invalid_true>0
  % keyboard
end
