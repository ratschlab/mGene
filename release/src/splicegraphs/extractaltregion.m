function [splicegraph, starts, stops] = extractaltregion(gene, region)

  start = region(1) ;
  stop  = region(2) ;
  
  exons = gene.splicegraph{1} ;
  introns= triu(gene.splicegraph{2}) ;

  idx=find(exons(1,:)<=stop & exons(2,:)>=start) ;
  
  %idx2 = [] ;
  %for i=1:length(idx)
  %  idx2 = [idx2 find(introns(:,idx(i)))' find(introns(idx(i),:))] ;
  %end ;
  %idx=unique([idx idx2]) ;
  
  splicegraph{1}=gene.splicegraph{1}(:,idx) ;
  splicegraph{2}=gene.splicegraph{2}(idx,idx) ;

  start_exons = idx(find(sum(introns(idx,idx))==0)) ;
  starts = exons(1,start_exons(1)):exons(2,start_exons(1));
  for i=2:length(start_exons)
    starts = intersect(starts,...
                       exons(1,start_exons(i)):exons(2,start_exons(i))) ;
  end ;
  if isempty(starts)
    starts = exons(1,idx(find(sum(introns(idx,idx))==0))) ;
  else
    starts = starts(1) ;
    splicegraph{1}(1,splicegraph{1}(1,:)<starts)=starts ;
  end ;
  
  stop_exons = idx(find(sum(introns(idx,idx)')==0)) ;
  stops = exons(1,stop_exons(1)):exons(2,stop_exons(1));
  for i=2:length(stop_exons)
    stops = intersect(stops,...
                      exons(1,stop_exons(i)):exons(2,stop_exons(i))) ;
  end ;
  if isempty(stops)
    stops = exons(2,idx(find(sum(introns(idx,idx)')==0))) ;
  else
    stops = stops(end) ;
    splicegraph{1}(2,splicegraph{1}(2,:)>stops)=stops ;
  end ;
  
  % merge exons, if they are identical
  [tmp,idx1,idx2]=unique(splicegraph{1}','rows') ;
  if size(tmp,1)~=size(splicegraph{1},2)
    for i=1:size(tmp,1)
      idx=find(idx2==i) ;
      if length(idx)>1,
        splicegraph{2}(idx,:)=repmat(sign(sum(splicegraph{2}(idx,:),1)),length(idx),1) ;
        splicegraph{2}(:,idx)=repmat(sign(sum(splicegraph{2}(:,idx),2)),1,length(idx)) ;
      end ;
    end ;
    splicegraph{1}=splicegraph{1}(:,idx1) ;
    splicegraph{2}=splicegraph{2}(idx1,idx1) ;
  end ;    
  