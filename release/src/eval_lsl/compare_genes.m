function [res_genes,map] = compare_genes(res_genes,annos,preds)

view_on = 0;
%map genes together
map=zeros(2,length(annos));

if isempty(annos)
  return;
end

num = 0;
rm_idx = [];
for i=1:length(preds)
  if isempty(preds(i).start)|isempty(preds(i).stop)
    rm_idx = [rm_idx i]
  end
end
preds(rm_idx) =[];
[p_starts,idx] = sort([preds.start]);
p_stops = [preds(idx).stop];
%p_ids = [preds(idx).id];
[a_starts,idx] = sort([annos.start]);
a_stops = [annos(idx).stop];
a_ids = [annos(idx).id];
for i=1:length(annos)
  idx = find(p_starts<=a_stops(i)&p_stops>=a_starts(i));
  if ~isempty(idx)
    idx2=num+[1:length(idx)] ;
    num=num+length(idx);
    map(2,idx2) = idx; %p_ids(idx);
    map(1,idx2) = a_ids(i);      
  end
end
map(:,num+1:end)=[];

wrong_ids = setdiff([preds.id],unique(map(2,:)));
missing_ids = setdiff([annos.id],unique(map(1,:)));

res_genes.num_obsv    = res_genes.num_obsv+length(annos);
res_genes.num_pred    = res_genes.num_pred+length(preds);
res_genes.wrong       = res_genes.wrong + (length(preds) - length(unique(map(2,:))));
res_genes.wrong_ids   = [res_genes.wrong_ids wrong_ids];

res_genes.missing     = res_genes.missing+ (length(annos) - length(unique(map(1,:))));
res_genes.missing_ids = [res_genes.missing_ids missing_ids];

status=inf;
for i=1:length(missing_ids)
  idx = find([annos.id]==missing_ids(i)) ;
  if ~isempty(idx) && isfield(annos, 'transcript_status') && ~isempty(annos(idx).transcript_status)
    status(i) = max(annos(idx).transcript_status);
  else
    status(i) = inf;
  end ;
end
res_genes.missing_conf = res_genes.missing_conf+sum(status==1);
res_genes.missing_partconf = res_genes.missing_partconf+sum(status==0);
res_genes.missing_unconf = res_genes.missing_unconf+sum(status==-1);


anno_ids = unique(map(1,:));
my_anno_ids = [annos.id] ;
num=0;
if isempty(anno_ids),
  anno_ids=[] ;
end ;

num_anno=[] ;
for i=anno_ids
  num=num+1;
  % fprintf('%i/%i\r',num, length(anno_ids))
  anno=annos(my_anno_ids==i);
  if ~(length(anno)==1),
    error('anno ids not unique') ;
  end ;
  num_anno(i) = sum(map(1,:)==i);
  corr=0;
  
  
  %% start loop over overlapping predicted genes
  for j=map(2,map(1,:)==i)
    pred = preds(j);
    
    for k=1 :length(anno.cds_exons)
      for l=1 :length(pred.cds_exons)
        % keyboard
        %if length(anno.transcript_status)>=k && anno.transcript_status(k)==1
          if isequal(pred.cds_exons{l}, anno.cds_exons{k})
            corr=1;
          end
        %elseif length(anno.transcript_status)>=k && (anno.transcript_status(k)==0|anno.transcript_status(k)==-1)%???
        %  if isequal(pred.cds_exons{l}, anno.cds_exons{k})
        %    % if all(ismember(anno.cds_exons{k},pred.cds_exons{l}))
        %    corr=1;
        %  end
        %end
      end %loop over pred transcripts
    end %loop over anno transcripts
    if view_on && ~corr
      anno.exons=pred.exons;
      viewsplicegraph(anno)
      pause
      % keyboard
    end
  end
  res_genes.num_corr=res_genes.num_corr+corr;
  %num_anno(i) = sum(map(1,:)==i);
end %%loop over anno_genes


res_genes.cut = res_genes.cut+ sum(num_anno>1); 
pred_ids = unique(map(2,:));
if ~isempty(pred_ids)
  for i=pred_ids
    num_pred(i) = sum(map(2,:)==i);
  end
  res_genes.merged = res_genes.merged + sum(num_pred(num_pred>1)); 
end



