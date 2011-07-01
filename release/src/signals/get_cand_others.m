function [POS,LABEL,CONF] = get_cand_others(P)
% [POS,LABEL,CONF] = get_cand_others(P)

debug = 0;  

genes = P.genes ;
seq = P.seq ;
Signal = P.Signal ;
info = P.info ;
strand = P.strand ;
signal_name = P.signal_name ;
%clear P

%% load signal consensus strings from file 
%% this should only apply for the polyA signal
if isempty(Signal.consensus) && isfield(Signal, 'consensus_fn') && ~isempty(Signal.consensus_fn)
  load(Signal.consensus_fn,'consensus')
  Signal.consensus = consensus;
end

%[Signal,info,info_names,sig] = initialize_cand_search(PAR,genes,strand,signal_name,debug);
[info_names,sig] = initialize_cand_search(Signal, info, genes, strand, signal_name, debug);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET ALL CONSENSUS SITES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(Signal.consensus)
  fprintf(1,'getting all consensus sites for: %s\n',Signal.consensus{:})
  POS =  find_consensus(seq,Signal,strand);
else
  fprintf(1,'no consensus given, all positions are potential candidates\n')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET ALL POSITIVE SITES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num = 0;
pos = zeros(1,length(seq));
conf = -ones(length(info_names),length(seq))*200;
gene_idx = zeros(1,length(seq));
for g_idx = 1:length(genes)
  gene = genes(g_idx) ;
  if ~isempty(gene.(signal_name))
    unique_pos = unique([gene.(signal_name){:}]);
    unique_pos2 = setdiff(unique_pos, pos(1:num)) ;
    if ~isequal(unique_pos2, unique_pos),
      %warning('identified other genes with same signal positions. Positions ignored')
      unique_pos = unique_pos2 ;
    end ;
    pos(num+[1:length(unique_pos)]) = unique_pos;
    gene_idx(num+[1:length(unique_pos)]) = repmat(g_idx,length(unique_pos),1);
    % === assertion fails for A. thaliana TSS !!! ===
    assert(length(unique(pos(1:num+length(unique_pos))))==length(pos(1:num+length(unique_pos))))
    conf_temp = zeros(length(info_names),length(unique_pos));
    for i=1:length(info_names)
      idx = find(gene.([signal_name '_info'])== ...
                 PAR.info_genes.(signal_name).(info_names{i}));
      pos_temp = [gene.(signal_name){idx}];
      [temp,idx1,idx2] = intersect(unique_pos,pos_temp);
      assert(length(idx) == length(pos_temp))
      conf_temp(i,idx1) = gene.([signal_name '_conf'])(idx(idx2));
    end
    conf(:,num+[1:length(unique_pos)]) = conf_temp; 
    num = num+length(unique_pos);
  end
end; %%%loop over genes
fprintf( 'found %d positive %s signals\n', num, signal_name );
pos(num+1:end) = [];
conf(:,num+1:end) = [];

if isempty(Signal.consensus)
  invalid_true = 0;
  neg_pos = find(~ismember([1:length(seq)],pos));
else
  invalid_true = length(setdiff(pos,POS));
  if invalid_true>10
    warning('%i positions did not exhibit the defined consensus', invalid_true) ;
  end
  pos = intersect(pos,POS); 
  neg_pos = setdiff(POS,pos);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLE DOWN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randn('state',24593854);
rand('state',24593854);

%neg_pos = sort(neg_pos) ;
fprintf(1,'Sampling down negatives:%d\n',Signal.sampling.downsample_neg);
if Signal.sampling.downsample_neg~=1
  var_ = 1/Signal.sampling.downsample_neg*1/5;
  %mean_ = Signal.sampling.mindist+1/Signal.sampling.downsample_neg;
  mean_ = 1/Signal.sampling.downsample_neg;
  num_new_pos = ceil((length(neg_pos)/mean_)*1.1);
  dist = abs(ceil(randn(num_new_pos,1).*var_+ mean_));
  %dist(dist<Signal.sampling.mindist)=dist(dist<Signal.sampling.mindist)+Signal.sampling.mindist;
  dist(dist==0) = deal(Signal.sampling.mindist);
  dd = cumsum(dist);
  neg_pos = intersect(neg_pos,dd);
  %keep = ceil(rand(1,length(neg_pos)*Signal.sampling.downsample_neg)*length(neg_pos));
  %keep = unique(keep);
  %neg_pos = neg_pos(keep);
elseif Signal.sampling.downsample_neg==1 & Signal.sampling.mindist>1
  remove = 1;
  while remove,
    dist = neg_pos(2:end)-neg_pos(1:end-1) ;
    idx = find(dist<Signal.sampling.mindist)+1;
    neg_pos(idx) = [];
    if isempty(idx)
      remove = 0;
    end
  end
end

POS = [pos neg_pos];
LABEL = [ones(1,length(pos)),-ones(1,length(neg_pos))]; 
COV1 = [conf zeros(size(conf,1),length(neg_pos))]; 


[temp,idx] = unique(POS);
assert(length(temp)==length(POS));
POS = POS(idx);
LABEL = LABEL(idx);
COV1 = COV1(:,idx);
if ~isequal(determine_engine, 'octave'), 
  assert(issorted(POS));
end ;


ALL = initialize_allconfs(POS,info_names);
ALL.COV1 = COV1;

  
fprintf('processing: ') ;
num_invalid_genes = 0;
for g_idx = 1:length(genes)
  if mod(g_idx, 100)==0,
    fprintf('%i/%i .. ', num, length(genes));
  end ;
  num = num+1;
  gene = genes(g_idx) ;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % GET BASIC STUFF, 
  % if candidate pos is in gene, get gene id, if gene is valid, if
  % candidate pos is in alternative region 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [ALL,is_invalid] = basic_cand_stuff(gene,ALL,sig,debug,signal_name);
  num_invalid_genes  = num_invalid_genes + is_invalid;
end %%%loop over genes
fprintf('Done\n') ;


if nargout>=3,
  CONF = set_CONF(Signal, ALL, sig, info_names, signal_name);
end ;


% assert(all(sum(Conf(1:3,LABEL==1),1)))
% assert(size(Conf,1)==length(Signal.Conf_names));
fprintf(1,'number of invalid genes: \t\t\t%i\n',num_invalid_genes)
fprintf(1,'number of true %s sites: \t\t\t%i\n',signal_name,sum(LABEL==1))
fprintf(1,'number of invalid true sites (no consensus): \t%i\n',invalid_true)
fprintf(1,'number of decoy sites:\t\t\t\t%i\n',length(LABEL))
% fprintf(1,'mean confidence:\t\t\t\t%.2f\n',mean(sum(Conf,1)));
