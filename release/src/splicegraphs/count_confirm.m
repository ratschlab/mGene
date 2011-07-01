function genes = count_confirm(genes)
% function genes = count_confirm(genes)
%
% Count the number of confirmations of alternative splicing events,
% and add to the genes structure.
%
% Written by: Cheng Soon Ong, 14 May 2007.




idx_alt=[] ;
for i=1:length(genes)
  if genes(i).is_alt_spliced
    idx_alt=[idx_alt, i] ;
  end
end

[idx_exon_skips, exon_exon_skips] = detect_exonskips(genes, idx_alt);
[idx_alt_5prime,exon_alt_5prime, idx_alt_3prime,exon_alt_3prime] = ...
    detect_altprime_pair(genes, idx_alt);
[idx_intron_reten,intron_intron_reten] = detect_intronreten(genes, idx_alt) ;

for ix = 1:length(genes)
  if mod(ix,1000)==0
    fprintf('.');
  end

  if ismember(ix,idx_exon_skips)
    tix = find(idx_exon_skips==ix);
    genes(ix).alt_exon_details.conf_count = zeros(2,0);
    for ixt = tix
      [conf_count0,conf_count1] = count_confirm_alt_exon(genes(ix),...
					exon_exon_skips(1,ixt),exon_exon_skips(2,ixt),exon_exon_skips(3,ixt));
      genes(ix).alt_exon_details.conf_count(1,end+1) = conf_count0;
      genes(ix).alt_exon_details.conf_count(2,end) = conf_count1;
    end
  end

  if ismember(ix,idx_intron_reten)
    tix = find(idx_intron_reten==ix);
    genes(ix).alt_intron_details.conf_count = zeros(2,0);
    for ixt = tix
      [conf_count0,conf_count1] = count_confirm_alt_intron(genes(ix),...
						  intron_intron_reten(1,ixt),intron_intron_reten(2,ixt));
      genes(ix).alt_intron_details.conf_count(1,end+1) = conf_count0;
      genes(ix).alt_intron_details.conf_count(2,end) = conf_count1;
    end
  end

  if ismember(ix,idx_alt_5prime)
    tix = find(idx_alt_5prime==ix);
    genes(ix).alt_5prime_details.conf_count = zeros(2,0);
    for ixt = tix
      if genes(ix).strands(1) == '+'
	idx_exon_l = exon_alt_5prime(ixt).fiveprimesites(1);
	idx_exon_r = exon_alt_5prime(ixt).threeprimesite;
	idx_exon_alt = exon_alt_5prime(ixt).fiveprimesites(2);
      else
	idx_exon_r = exon_alt_5prime(ixt).fiveprimesites(1);
	idx_exon_l = exon_alt_5prime(ixt).threeprimesite;
	idx_exon_alt = exon_alt_5prime(ixt).fiveprimesites(2);
      end
      [conf_count0,conf_count1] = count_confirm_alt_5prime(genes(ix),idx_exon_l,idx_exon_r,idx_exon_alt);
      genes(ix).alt_5prime_details.conf_count(1,end+1) = conf_count0;
      genes(ix).alt_5prime_details.conf_count(2,end) = conf_count1;
    end
  end

  if ismember(ix,idx_alt_3prime)
    tix = find(idx_alt_3prime==ix);
    genes(ix).alt_3prime_details.conf_count = zeros(2,0);
    for ixt = tix
      if genes(ix).strands(1) == '+'
	idx_exon_l = exon_alt_3prime(ixt).fiveprimesite;
	idx_exon_r = exon_alt_3prime(ixt).threeprimesites(1);
	idx_exon_alt = exon_alt_3prime(ixt).threeprimesites(2);
      else
	idx_exon_r = exon_alt_3prime(ixt).fiveprimesite;
	idx_exon_l = exon_alt_3prime(ixt).threeprimesites(1);
	idx_exon_alt = exon_alt_3prime(ixt).threeprimesites(2);
      end
      [conf_count0,conf_count1] = count_confirm_alt_3prime(genes(ix),idx_exon_l,idx_exon_r,idx_exon_alt);
      genes(ix).alt_3prime_details.conf_count(1,end+1) = conf_count0;
      genes(ix).alt_3prime_details.conf_count(2,end) = conf_count1;
    end
  end

end
