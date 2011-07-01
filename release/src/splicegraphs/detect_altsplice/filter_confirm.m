function [eec,eic,e5c,e3c] = filter_confirm(genes,ee,ei,e5,e3,min_confirmed_pos0,min_confirmed_pos1)
% function [eec,eic,e5c,e3c] = filter_confirm(genes,ee,ei,e5,e3,min_confirmed_pos0,min_confirmed_pos1)
%
% filter list of genes to only have events with large numbers of confirmed alternative splicing

if nargin < 6
  min_confirmed_pos0 = 3;
  min_confirmed_pos1 = 3;
end

% find the alternatively spliced genes
idx_alt=[] ;
for i=1:length(genes)
  if genes(i).is_alt_spliced
    idx_alt=[idx_alt, i] ;
  end
end

[idx_exon_skips, exon_exon_skips] = detect_exonskips(genes, idx_alt);
[idx_alt_5prime,exon_alt_5prime, idx_alt_3prime,exon_alt_3prime] = ...
    detect_altprime(genes, idx_alt);
[idx_intron_reten,intron_intron_reten] = detect_intronreten(genes, idx_alt) ;

eec = [];
for ix = ee
  cur_idx = find(ix==idx_exon_skips);
  for idx = cur_idx
    [conf_count0,conf_count1]=count_confirm_alt_exon(genes(ix),exon_exon_skips(1,idx),...
						     exon_exon_skips(2,idx),exon_exon_skips(3,idx));
    if (conf_count0>=min_confirmed_pos0) && (conf_count1>=min_confirmed_pos1)
      eec(end+1) = ix;
    end
  end
end

eic = [];
for ix = ei
  cur_idx = find(ix==idx_intron_reten);
  for idx = cur_idx
    [conf_count0,conf_count1] = count_confirm_alt_intron(genes(ix),intron_intron_reten(1,idx),...
						  intron_intron_reten(2,idx));
    if (conf_count0>=min_confirmed_pos0) && (conf_count1>=min_confirmed_pos1)
      eic(end+1) = ix;
    end
  end
end

e3c= [];
for ix = e3
  cur_idx = find(ix==idx_alt_3prime);
  for idx = cur_idx
    if genes(ix).strands(1) == '+'
      [conf_count0,conf_count1] = count_confirm_alt_3prime(genes(ix),exon_alt_3prime(idx).fiveprimesite,...
						  exon_alt_3prime(idx).threeprimesites(1),...
						  exon_alt_3prime(idx).threeprimesites(end));
      if (conf_count0>=min_confirmed_pos0) && (conf_count1>=min_confirmed_pos1)
	e3c(end+1) = ix;
	continue;
      end
      for ixs = 2:length(exon_alt_3prime(idx).threeprimesites)
	[conf_count0,conf_count1] = count_confirm_alt_3prime(genes(ix),exon_alt_3prime(idx).fiveprimesite,...
						  exon_alt_3prime(idx).threeprimesites(ixs),...
						  exon_alt_3prime(idx).threeprimesites(ixs-1));
	if (conf_count0>=min_confirmed_pos0) && (conf_count1>=min_confirmed_pos1)
	  e3c(end+1) = ix;
	  continue;
	end
      end
    else
      [conf_count0,conf_count1] = count_confirm_alt_3prime(genes(ix),exon_alt_3prime(idx).threeprimesites(1),...
						  exon_alt_3prime(idx).fiveprimesite,...
						  exon_alt_3prime(idx).threeprimesites(end));
      if (conf_count0>=min_confirmed_pos0) && (conf_count1>=min_confirmed_pos1)
	e3c(end+1) = ix;
	continue;
      end
      for ixs = 2:length(exon_alt_3prime(idx).threeprimesites)
	[conf_count0,conf_count1] = count_confirm_alt_3prime(genes(ix),exon_alt_3prime(idx).threeprimesites(ixs),...
						  exon_alt_3prime(idx).fiveprimesite,...
						  exon_alt_3prime(idx).threeprimesites(ixs-1));
	if (conf_count0>=min_confirmed_pos0) && (conf_count1>=min_confirmed_pos1)
	  e3c(end+1) = ix;
	  continue;
	end
      end
    end % if genes(ix).strands(1) == '+'
  end
end



e5c= [];
for ix = e5
  cur_idx = find(ix==idx_alt_5prime);
  for idx = cur_idx
    if genes(ix).strands(1) == '+'
      [conf_count0,conf_count1] = count_confirm_alt_5prime(genes(ix),exon_alt_5prime(idx).fiveprimesites(1),...
						  exon_alt_5prime(idx).threeprimesite,...
						  exon_alt_5prime(idx).fiveprimesites(end));
      if (conf_count0>=min_confirmed_pos0) && (conf_count1>=min_confirmed_pos1)
	e5c(end+1) = ix;
	continue;
      end
      for ixs = 2:length(exon_alt_5prime(idx).threeprimesites)
	[conf_count0,conf_count1] = count_confirm_alt_5prime(genes(ix),exon_alt_5prime(idx).fiveprimesites(ixs),...
						  exon_alt_5prime(idx).threeprimesite,...
						  exon_alt_5prime(idx).fiveprimesites(ixs-1));
	if (conf_count0>=min_confirmed_pos0) && (conf_count1>=min_confirmed_pos1)
	  e5c(end+1) = ix;
	  continue;
	end
      end
    else
      [conf_count0,conf_count1] = count_confirm_alt_5prime(genes(ix),exon_alt_5prime(idx).threeprimesite,...
						  exon_alt_5prime(idx).fiveprimesites(1),...
						  exon_alt_5prime(idx).fiveprimesites(end));
      if (conf_count0>=min_confirmed_pos0) && (conf_count1>=min_confirmed_pos1)
	e5c(end+1) = ix;
	continue;
      end
      for ixs = 2:length(exon_alt_5prime(idx).threeprimesites)
	[conf_count0,conf_count1] = count_confirm_alt_5prime(genes(ix),exon_alt_5prime(idx).threeprimesite,...
						  exon_alt_5prime(idx).fiveprimesites(ixs),...
						  exon_alt_5prime(idx).fiveprimesites(ixs-1));
	if (conf_count0>=min_confirmed_pos0) && (conf_count1>=min_confirmed_pos1)
	  e5c(end+1) = ix;
	  continue;
	end
      end
    end % if genes(ix).strands(1) == '+'
  end
end

