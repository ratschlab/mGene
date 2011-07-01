
if 0,
  level_e = zeros(1,length(genes));
  level_i = zeros(1,length(genes));
  for ix = 1:length(genes)
    if (mod(ix,100)==0)
      fprintf(1,'.');
    end
    [level_e(ix), level_i(ix)] = detectsplicegraph(genes(ix));
  end
  
  fprintf(1,'\n\nTotal alternatively spliced:\t\t\t%d\n',...
	  sum(or(level_i>1,level_e>1)));
  fprintf(1,'Total constitutively spliced:\t\t\t%d\n',...
	  sum(and(level_i<=1,level_e<=1)));
  
  [dum,idx_alt] = find(or(level_i>1,level_e>1));
  [dum,idx_con] = find(and(level_i<=1,level_e<=1));
  
  
  
  disp('Sorting exons...');
  for ix=idx_alt
    [dummy,exon_order] = sort(genes(ix).splicegraph{1}(1,:),2,'ascend');
    genes(ix).splicegraph{1} = genes(ix).splicegraph{1}(:,exon_order);
    genes(ix).splicegraph{2} = genes(ix).splicegraph{2}(exon_order,exon_order);
  end
end;








idx_multiple_skips = [];
exon_multiple_skips = [];
%for ix=idx_alt
for ix = 1:length(genes)
  fprintf('%i\r', ix) ;

  num_exons = size(genes(ix).splicegraph{1},2);
  edges = genes(ix).splicegraph{2};
  labels = repmat([1:num_exons]',1,num_exons);
  for exon_idx_first = 1:num_exons-3
    for exon_idx_last = (exon_idx_first+3):num_exons
      if (edges(exon_idx_first,exon_idx_last))
	
	
	
	% find all pairs shortest path
	exist_path = triu(edges);
	exist_path(exon_idx_first,exon_idx_last)=0;
	exist_path(exist_path==0)=inf ;
	for i=1:num_exons; exist_path(i,i)=0 ; end ;

	long_exist_path = -triu(edges);
	long_exist_path(exon_idx_first,exon_idx_last)=0;
	long_exist_path(long_exist_path==0)=inf ;
	for i=1:num_exons; long_exist_path(i,i)=0 ; end ;

	
	path = ~isinf(exist_path).*labels;
	long_path = ~isinf(long_exist_path).*labels;
	
	for k=1:num_exons,
	  for i=1:num_exons,
	    idx=find(exist_path(i,k)+exist_path(k,:)<exist_path(i,:));
	    exist_path(i,idx)=exist_path(i,k)+exist_path(k,idx) ;
	    path(i,idx) = path(k,idx) ;

	    idx=find(long_exist_path(i,k)+long_exist_path(k,:)<long_exist_path(i,:));
	    long_exist_path(i,idx)=long_exist_path(i,k)+long_exist_path(k,idx);
	    long_path(i,idx) = long_path(k,idx);
	  end 
	end 

	temp_ix = ~isinf(long_exist_path);
	long_exist_path(temp_ix) = -long_exist_path(temp_ix);
	
	if exist_path(exon_idx_first,exon_idx_last)>2 & ~isinf(exist_path(exon_idx_first,exon_idx_last))
	  backtrace = path(exon_idx_first,exon_idx_last);
	  while (backtrace(end) > exon_idx_first)
	    backtrace = [backtrace,path(exon_idx_first,backtrace(end))];
	  end
	  backtrace = backtrace(1:end-1);
	  backtrace = backtrace(length(backtrace):-1:1);
	  idx_multiple_skips = [idx_multiple_skips,repmat(ix,1,length(backtrace))];
	  exon_multiple_skips = [exon_multiple_skips,backtrace];
	elseif long_exist_path(exon_idx_first,exon_idx_last)>2 & ~isinf(long_exist_path(exon_idx_first,exon_idx_last))
	  backtrace = long_path(exon_idx_first,exon_idx_last);
	  while (backtrace(end) > exon_idx_first)
	    backtrace = [backtrace,long_path(exon_idx_first,backtrace(end))];
	  end
	  backtrace = backtrace(1:end-1);
	  backtrace = backtrace(length(backtrace):-1:1);
	  idx_multiple_skips = [idx_multiple_skips,repmat(ix,1,length(backtrace))];
	  exon_multiple_skips = [exon_multiple_skips,backtrace];
	end
	%if exist_path(exon_idx_first,exon_idx_last) < ...
	%      long_exist_path(exon_idx_first,exon_idx_last)
	%  disp(ix);
	%end
	%if exist_path(exon_idx_first,exon_idx_last) > ...
	%      long_exist_path(exon_idx_first,exon_idx_last)
	%  keyboard
	%end
	
      end
    end
  end
end
fprintf(1,'Number of multiple exon skips:\t\t\t\t%d\n',...
	length(idx_multiple_skips));




%for ix=1:length(idx_multiple_skips);
%  disp(idx_multiple_skips(ix));
%  viewsplicegraph_compare(genes(idx_multiple_skips(ix)),genes(idx_multiple_skips(ix)).splicegraph{1}(:,exon_multiple_skips(ix)));
%  pause;
%end;
