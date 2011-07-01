function [level,cds_level]  = cnt_obsv_pred_corr(level,cds_level,anno,pred,strand,dist)

if ~isfield(anno,'all_gene_ids')
  anno.all_gene_ids =[];
  anno.cds_gene_ids =[];
end
if ~isfield(pred,'all_gene_ids')
  pred.all_gene_ids =[];
  pred.cds_gene_ids =[];
end
if ~isfield(anno,'status_all')
  anno.status_all =[];
  anno.status_cds =[];
end

  
if ~isempty(cds_level)
  level     = cnt_tmp(level,anno.all,pred.all,anno.all_gene_ids,pred.all_gene_ids,anno.status_all);
  cds_level = cnt_tmp(cds_level,anno.cds,pred.cds,anno.cds_gene_ids,pred.cds_gene_ids,anno.status_cds);
else
  names = {'tis','cdsStop','acc','don'};
  for i=1:length(names)
    if isfield(anno,names{i})
      names{i};
      level.(names{i})     = cnt_tmp(level.(names{i}),anno.(names{i}),pred.(names{i}),anno.all_gene_ids,pred.all_gene_ids,anno.status_all);
    end
  end
  names = {'tss','cleave'};
  for i=1:length(names)
    if isfield(anno,names{i})
      names{i};
      level.(names{i})     = cnt_tmp(level.(names{i}),anno.(names{i}),pred.(names{i}),anno.all_gene_ids,pred.all_gene_ids,anno.status_all,strand,dist);
    end
  end
end


function level = cnt_tmp(level,anno,pred,anno_gene_ids,pred_gene_ids,anno_status,strand,dist)
  

level.num_obsv = level.num_obsv + size(anno,1);
level.num_pred = level.num_pred + size(pred,1);

if min(size(anno))==1 && min(size(pred))<=1 &&  ~iscell(anno) % nucleotides % signals
  level.num_corr = level.num_corr + length(intersect(anno,pred));  
  
  if exist('dist')
    D=[];I=[];
    for p=pred'
      [tmp,nb] = min(abs(anno-p));
      if strand=='+'
        D = [D anno(nb)-p];
      else
        D = [D -anno(nb)+p];
      end
      I=[I nb];
    end
    level.dist=[level.dist D];
    for d=dist
      idx = abs(D)<=d;
      num = length(unique(I(idx)));
      level.(['num_corr' num2str(d)]) = level.(['num_corr' num2str(d)])+num;
    end
  end
elseif ~iscell(anno)  % exons
  level.num_corr = level.num_corr + size(intersect(anno,pred,'rows'),1);  
  for i=1:size(anno,1)
    if isempty(pred)||~any(pred(:,1)<=anno(i,2)&pred(:,2)>=anno(i,1))
      level.ME = level.ME +1 ;
      level.ME_gene_id = [level.ME_gene_id anno_gene_ids(i)];
    end
  end
  for i=1:size(pred,1)
    if isempty(anno), continue; end
    if ~any(anno(:,1)<=pred(i,2)&anno(:,2)>=pred(i,1))
      level.WE = level.WE +1 ;
      level.WE_gene_id = [level.WE_gene_id pred_gene_ids(i)];
    end
  end
  
else %transcripts
  corr = 0;
  if length(anno)==0
	return
  end
  a_starts=[] ;
  a_stops=[] ;
  a_idx=[] ;
  for j=1:length(anno)
    if ~isempty(anno{j}),
      a_starts(end+1) = min(anno{j}(:));
      a_stops(end+1) = max(anno{j}(:));
      a_idx(end+1) = j ;
    end ;
  end
  for j=1:length(pred)
    start = min(pred{j}(:));
    stop = max(pred{j}(:));
    if isempty(start) || isempty(stop)
      continue
    end
    idx = find(a_starts==start & a_stops==stop) ;
    if isempty(idx)
      continue
    end
    
    for i=idx       
      if isequal(anno{a_idx(i)}, pred{j})
        corr=corr+1;
        break;
      end
    end
  end
  level.num_corr = level.num_corr+corr;
end


