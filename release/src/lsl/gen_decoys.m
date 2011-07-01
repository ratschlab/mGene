function decoys = gen_decoys(block, model, max_num)
% decoys = gen_decoys(block, model, max_num)

 
gene_idx = setdiff(unique(block.truth(1).segments(:,4)),0);
for j=1:length(gene_idx)
  idx = find(block.truth(1).segments(:,4)==gene_idx(j));
  genes{j} = block.truth(1).segments(idx,1:4) ;
end
decoys.genes = genes ;
%decoys.genes = block.truth(1).genes ;

modified = 0 ;
for i = 1:length(decoys.genes)
  genes = decoys.genes{i} ;
  idx = find(genes(:,3)==model.segments.cds_exon) ;
  idx = idx(2:end-1) ;
  idx1 = find(mod(genes(idx,2)-genes(idx,1),3)==0);
  if length(idx1)>0,
    decoys.genes{i}(idx(idx1)-1,2) = decoys.genes{i}(idx(idx1)+1,2) ;
    decoys.genes{i}(idx(idx1):idx(idx1)+1,:) = [] ;
    modified = 1 ;
  end ;
end ;

if ~modified, 
  decoys = [] ; 
  return ; 
end ;

idx = find(block.truth(1).segments(:,3)==model.segments.intergenic);
if isfield(model.segments, 'intercistronic')
  idx = [idx find(block.truth(1).segments(:,3)==model.segments.intercistronic)];
end
decoys.segments = [] ;
for i = 1:length(decoys.genes)
  decoys.segments = [decoys.segments; block.truth(1).segments(idx(i),1:4)];
  decoys.segments = [decoys.segments; decoys.genes{i}] ;
end ;
decoys.segments = [decoys.segments;block.truth(1).segments(end,1:4)];




[decoys.path, decoys.pos_idx, decoys.pos] = segmentation2path(decoys, block, model);

decoys = orderfields(decoys) ;
return ;




















truepos = block.truth.exons{1}(:,1:2) ;
truedescr = block.truth.exons{1}(:,3) ;
content_type =  unique(truedescr) ;

paths = {} ;
descr = {} ;
num = 0 ;

% simulate gene skip
gene_id = find(truedescr == getfield(model.contents,'intergenic')) ;
true_intergenic = truepos(idx) ;
for i = 1:size(gene_id,1)-2,
  num = num+1 ;
  paths{num} = [truepos(1:gene_id(i)-1,:);
                truepos(gene_id(i),1) truepos(gene_id(i)+1,2) 
                truepos(gene_id(i)+2:end,:)] ;
  descr{num} = [truedescr(1:gene_id(i));
                truedescr(gene_id(i)+2:end)] ;
end ;


% simulate merged gene 
for i = 1:size(gene_id,1)-1,
  num = num+1 ;
  paths{num} = [truepos(1:gene_id(i)-1,:);
                truepos(gene_id(i+1):end,:)] ;
  descr{num} = [truedescr(1:gene_id(i)-1);truedescr(gene_id(i+1):end)] ;
end ;


% simulate intron skip
for g = 1:length(gene_id)
  truepos_g = truepos(gene_id(g):gene_id(g+1)-1,:);
  for i = 1:size(truepos_g,1)-2,
    if truedescr(i) == truedescr(i+1)
      num = num+1 ;
      paths{num} = truepos ;
      
      paths{num} = [truepos(1:i-1,:);
                    truepos(i,1) truepos(i+1,2) 
                    truepos(i+2:end,:)] ;
      descr{num} = [truedescr(1:i);
                    truedescr(i+2:end)] ;
    end
  end ;
end

% simulate exon skip
for i=2:size(truepos,1)-1,
  num=num+1 ;
  paths{num}=[truepos(1:i-1,:);
              truepos(i+1:end,:)] ;
  descr{num} = [truedescr(1:i-1);truedescr(i+1:end)] ;
end ;

% simulate alternative donor
for w = 1:length(unique(content_type))
  truepos_w = truepos(:,truedescr==w)
  for i = 1:size(truepos_w,1)-1,
    num = num+1 ;
    paths{num} = truepos ;
    descr{num} = truedescr
    range = truepos_w(i,1)+3:truepos_w(i+1,1) ;
    idx = range(strfind(seq(range),'gt')) ;
    idx = idx(idx~=truepos(i,2)) ;
    idx2 = randperm(length(idx)) ;
    if isempty(idx2),
      num = num-1 ;
    else
      paths{num}(i,2) = idx(idx2(1)) ;
    end ;
  end ;
end

% simulate alternative acceptor
for w = 1:length(unique(content_type))
  truepos_w = truepos(:,truedescr==w)
  for i = 1:size(truepos,1)-1,
    num = num+1 ;
    paths{num} = truepos ;
    descr{num} = truedescr
    range = truepos_w(i,2):truepos_w(i+1,2)-3 ;
    idx = range(strfind(seq(range),'ag'))+2 ;
    idx = idx(idx~=truepos(i+1,1)) ;
    idx2 = randperm(length(idx)) ;
    if isempty(idx2),
      num = num-1 ;
    else
      paths{num}(i+1,1)=idx(idx2(1)) ;
    end ;
  end ;
end

% simulate alternative starts 
for w = 1:length(unique(content_type))
  truepos_w = truepos(:,truedescr==w)
  for i = 1:size(truepos_w,1)-1,
    num = num+1 ;
    paths{num} = truepos ;
    descr{num} = truedescr
    range = truepos_w(i,1)+3:truepos_w(i+1,1) ;
    idx = range(strfind(seq(range),'gt')) ;
    idx = idx(idx~=truepos(i,2)) ;
    idx2 = randperm(length(idx)) ;
    if isempty(idx2),
      num = num-1 ;
    else
      paths{num}(i,2) = idx(idx2(1)) ;
    end ;
  end ;
end


% simulate alternative PolyA

% simulate alternative TIS

% simulate alternative STOP


for i=1:num
  assert(all(paths{i}(:,1)<paths{i}(:,2))) ;
  assert(all(paths{i}(1:end-1,2)<paths{i}(2:end,1))) ;
end ;

idx = randperm(num) ;
paths = paths(idx(1:min(max_num,num))) ;


for i=1:length(paths)
  decoys(i).exons = paths{i} ;
  decoys(i).path = [] ;
end ;
