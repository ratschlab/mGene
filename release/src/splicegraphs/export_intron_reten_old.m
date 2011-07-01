% Generate a matlab struct containing the information required for
% classification of intron retentions.

%load genes.all.graph.mat


%countspliceform

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
  [dummy,exon_order] = sort(genes(ix).splicegraph{1}(1,:));
  genes(ix).splicegraph{1} = genes(ix).splicegraph{1}(:,exon_order);
  genes(ix).splicegraph{2} = genes(ix).splicegraph{2}(exon_order,exon_order);
end



genes_intron_reten = [];
% construct the output
genes_intron_reten(1).id = 0;
genes_intron_reten(1).name = 'a';
genes_intron_reten(1).chr = 'a';
genes_intron_reten(1).strand = 'a';
genes_intron_reten(1).label = 99;
genes_intron_reten(1).exon_5prime = [0 0];
genes_intron_reten(1).exon_3prime = [0 0];
genes_intron_reten(1).exon_intron = [0 0];
genes_intron_reten(1).confirm_by = 0;
	  



% Generate positive examples
idx_intron_reten = [];
for ix=idx_alt
  num_exons = size(genes(ix).splicegraph{1},2);
  vertices = genes(ix).splicegraph{1};
  edges = genes(ix).splicegraph{2};
  for exon_idx = 1:num_exons-2
    for exon_idx1 = (exon_idx+1):(num_exons-1)
      for exon_idx2 = (exon_idx1+1):num_exons
	if (edges(exon_idx,exon_idx2)&&...
	    (vertices(1,exon_idx1)<vertices(2,exon_idx))&&...
	    (vertices(2,exon_idx1)>vertices(1,exon_idx2))&&...
	    (vertices(2,exon_idx)<vertices(2,exon_idx1))&&...
	    (vertices(1,exon_idx1)<vertices(1,exon_idx2)))
	  idx_intron_reten = [idx_intron_reten,ix];
	  
	  conf_count = 0;
	  for idx1 = 1:length(genes(ix).exons)
	    for idx2 = 1:size(genes(ix).exons{idx1},1)
	      if (vertices(2,exon_idx) > genes(ix).exons{idx1}(idx2,1))&&...
		    (vertices(1,exon_idx2) < genes(ix).exons{idx1}(idx2,2))
		conf_count = conf_count + 1;
	      end
	    end
	  end
	  
	  
	  % construct the output
	  genes_intron_reten(end+1).id = ix;
	  genes_intron_reten(end).name = genes(ix).name;
	  genes_intron_reten(end).chr = genes(ix).chr;
	  genes_intron_reten(end).strand = genes(ix).strands(1);
	  genes_intron_reten(end).label = 1;
	  if (genes(ix).strands(1) == '+')
	    genes_intron_reten(end).exon_5prime = [vertices(1,exon_idx) vertices(2,exon_idx)];
	    genes_intron_reten(end).exon_3prime = [vertices(1,exon_idx2) vertices(2,exon_idx2)];
	  else
	    genes_intron_reten(end).exon_3prime = [vertices(1,exon_idx) vertices(2,exon_idx)];
	    genes_intron_reten(end).exon_5prime = [vertices(1,exon_idx2) vertices(2,exon_idx2)];
	  end
	  genes_intron_reten(end).exon_intron = [vertices(1,exon_idx1) vertices(2,exon_idx1)];
	  genes_intron_reten(end).confirm_by = conf_count;
	  
	  
	  
	elseif  (edges(exon_idx1,exon_idx2)&&...
	       (vertices(1,exon_idx)<vertices(2,exon_idx1))&&...
	       (vertices(2,exon_idx)>vertices(1,exon_idx2))&&...
	       (vertices(2,exon_idx)>vertices(2,exon_idx1))&&...
	       (vertices(1,exon_idx)<vertices(1,exon_idx2)))
	  idx_intron_reten = [idx_intron_reten,ix];

	  conf_count = 0;
	  for idx1 = 1:length(genes(ix).exons)
	    for idx2 = 1:size(genes(ix).exons{idx1},1)
	      if (vertices(2,exon_idx1) > genes(ix).exons{idx1}(idx2,1))&&...
		    (vertices(1,exon_idx2) < genes(ix).exons{idx1}(idx2,2))
		conf_count = conf_count + 1;
	      end
	    end
	  end
	  
	  
	  % construct the output
	  genes_intron_reten(end+1).id = ix;
	  genes_intron_reten(end).name = genes(ix).name;
	  genes_intron_reten(end).chr = genes(ix).chr;
	  genes_intron_reten(end).strand = genes(ix).strands(1);
	  genes_intron_reten(end).label = 1;
	  if (genes(ix).strands(1) == '+')
	    genes_intron_reten(end).exon_5prime = [vertices(1,exon_idx1) vertices(2,exon_idx1)];
	    genes_intron_reten(end).exon_3prime = [vertices(1,exon_idx2) vertices(2,exon_idx2)];
	  else
	    genes_intron_reten(end).exon_3prime = [vertices(1,exon_idx1) vertices(2,exon_idx1)];
	    genes_intron_reten(end).exon_5prime = [vertices(1,exon_idx2) vertices(2,exon_idx2)];
	  end	    
	  genes_intron_reten(end).exon_intron = [vertices(1,exon_idx) vertices(2,exon_idx)];
	  genes_intron_reten(end).confirm_by = conf_count;
	end
      end
    end
  end
end
fprintf(1,'Number of intron retentions:\t\t\t\t%d\n',...
	length(idx_intron_reten));



genes_intron_reten = genes_intron_reten(2:end);




  


% Also don't use introns which have exons incompletely in them
% Exon in intronic region, may be intron retention, but not complete.
idx_exon_intron = [];
for ix=idx_alt
  num_exons = size(genes(ix).splicegraph{1},2);
  vertices = genes(ix).splicegraph{1};
  edges = genes(ix).splicegraph{2};
  for exon_idx = 1:num_exons-2
    for exon_idx1 = (exon_idx+1):(num_exons-1)
      for exon_idx2 = (exon_idx1+1):num_exons
	if (edges(exon_idx1,exon_idx2)&&...
	    (sum(edges(exon_idx,exon_idx1:num_exons))==0)&&...
	    (vertices(2,exon_idx)>vertices(2,exon_idx1))&&...
	    (vertices(2,exon_idx)<=vertices(1,exon_idx2)))||...
	      (edges(exon_idx,exon_idx2)&&...
	       (sum(edges(exon_idx1,exon_idx2:num_exons))==0)&&...
	       (vertices(2,exon_idx1)>vertices(2,exon_idx))&&...
	       (vertices(2,exon_idx1)<=vertices(1,exon_idx2)))||...
	      (edges(exon_idx,exon_idx2)&&...
	       (sum(edges(1:exon_idx,exon_idx1))==0)&&...
	       (vertices(1,exon_idx1)>=vertices(2,exon_idx))&&...
	       (vertices(1,exon_idx1)<vertices(1,exon_idx2)))
	  idx_exon_intron = [idx_exon_intron,ix];
	end
      end
    end
  end
end
fprintf(1,'Number of incomplete exons in intronic regions:\t\t%d\n',...
	length(idx_exon_intron));



% Generate negative examples
idx_negative = setdiff([1:length(genes)],idx_intron_reten);
idx_negative = setdiff([1:length(genes)],idx_exon_intron);
fprintf(1,'Number of negative examples:\t\t\t\t%d\n',...
	length(idx_negative));

for ix=idx_negative
  num_exons = size(genes(ix).splicegraph{1},2);
  vertices = genes(ix).splicegraph{1};
  edges = genes(ix).splicegraph{2};

  for exon_idx1 = 1:(num_exons-1)
    for exon_idx2 = (exon_idx1+1):num_exons
      if (edges(exon_idx1,exon_idx2))
	
	
	
	conf_count = 0;
	for idx1 = 1:length(genes(ix).exons)
	  for idx2 = 1:size(genes(ix).exons{idx1},1)-1
	    if (vertices(2,exon_idx1) == genes(ix).exons{idx1}(idx2,2))&&...
		  (vertices(1,exon_idx2) == genes(ix).exons{idx1}(idx2+1,1))
	      conf_count = conf_count + 1;
	    end
	  end
	end
	  

	
	
	% construct the output
	genes_intron_reten(end+1).id = ix;
	genes_intron_reten(end).name = genes(ix).name;
	genes_intron_reten(end).chr = genes(ix).chr;
	genes_intron_reten(end).strand = genes(ix).strands(1);
	genes_intron_reten(end).label = -1;
	if (genes(ix).strands(1) == '+')
	  genes_intron_reten(end).exon_5prime = [vertices(1,exon_idx1) vertices(2,exon_idx1)];
	  genes_intron_reten(end).exon_3prime = [vertices(1,exon_idx2) vertices(2,exon_idx2)];
	else
	  genes_intron_reten(end).exon_3prime = [vertices(1,exon_idx1) vertices(2,exon_idx1)];
	  genes_intron_reten(end).exon_5prime = [vertices(1,exon_idx2) vertices(2,exon_idx2)];
	end
	genes_intron_reten(end).exon_intron = [];
	genes_intron_reten(end).confirm_by = conf_count;
      end
    end
  end
  
end
