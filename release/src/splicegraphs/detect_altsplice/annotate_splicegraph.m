function genes = annotate_splicegraph(genes)
% function annotate_splicegraph(info,origpath)
%
% annotate the genes structure with counts for:
% - exon skipping
% - intron retention
% - alternative donor (5 prime)
% - alternative acceptor (3 prime)
% - alternative transcription start
% - alternative transcription end
%
% Written by: Cheng Soon Ong, 14 May 2007
% Modified on: 31 May 2007 to include transcription starts and stops


%addpath(sprintf('%s/utils', origpath)) ;
%addpath(sprintf('%s/splicegraphs', origpath)) ;
%addpath(sprintf('%s/splicegraphs/detect_altsplice', origpath)) ;
%addpath(sprintf('%s/gff', origpath)) ;
%addpath(sprintf('%s/gff/read_sequences', origpath)) ;
%addpath /fml/ag-raetsch/share/software/matlab_tools/utils

%genome_info = init_genome(info) ;
%fprintf('directory: %s\n', genome_info.basedir) ;
%load(sprintf('%s/confirmed_sequences.mat', genome_info.basedir), 'genes') ;
%fprintf('loaded confirmed_sequences.mat\n\n') ;

if ~isfield(genes, 'strands'),
  for i=1:length(genes),
    genes(i).strands = repmat(genes(i).strand, 1, length(genes(i).transcripts)) ;
  end ;
end ;

if isfield(genes,'alt_3prime_details') || isfield(genes,'alt_5prime_details')
  warning('annotate_splicegraph adds to alt_3prime_details and alt_5prime_details,\nrunning more than one creates inconsistencies');
  return;
end

idx_alt=[] ;
for i=1:length(genes)
  if genes(i).is_alt_spliced
    idx_alt=[idx_alt, i] ;
  end
  if ~isfield(genes(i),'strands') || isempty(genes(i).strands)
    genes(i).strands(1) = genes(i).strand;
  end
end

[idx_exon_skips, exon_exon_skips] = detect_exonskips(genes, idx_alt);
[idx_alt_5prime,exon_alt_5prime, idx_alt_3prime,exon_alt_3prime] = ...
    detect_altprime_pair(genes, idx_alt);
[idx_intron_reten,intron_intron_reten] = detect_intronreten(genes, idx_alt) ;
[idx_xor_exons, exon_xor_exons] = detect_xorexons(genes, idx_alt) ;
[idx_alt_tstart, exon_alt_tstart] = detect_alttstart(genes) ;
idx_alt_tstart = [idx_alt_tstart{:}];
[idx_alt_tend, exon_alt_tend] = detect_alttend(genes) ;
idx_alt_tend = [idx_alt_tend{:}];

%keyboard
for ix = 1:length(genes)
    
  if mod(ix,1000)==0
    fprintf('.');
  end

  exons_is_alt = zeros(1,size(genes(ix).splicegraph{1},2));

  % exon skips
  if ismember(ix,idx_exon_skips)
    tix = find(idx_exon_skips==ix);
    genes(ix).alt_exon = length(tix);
    genes(ix).alt_exon_details.exons = exon_exon_skips(:,tix);
    exons_is_alt(exon_exon_skips(:,tix)) = exons_is_alt(exon_exon_skips(:,tix)) + 1;
  else
    genes(ix).alt_exon = 0;
  end

  % intron retentions
  if ismember(ix,idx_intron_reten)
    tix = find(idx_intron_reten==ix);
    genes(ix).alt_intron = length(tix);
    genes(ix).alt_intron_details.exons = intron_intron_reten(:,tix);
    exons_is_alt(intron_intron_reten(:,tix)) = exons_is_alt(intron_intron_reten(:,tix)) + 1;
  else
    genes(ix).alt_intron = 0;
  end

  % alternative 5prime
  if ismember(ix,idx_alt_5prime)
    genes(ix).alt_5prime_details.exons = zeros(3,0);
    tix = find(idx_alt_5prime==ix);
    genes(ix).alt_5prime = length(tix);
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
      genes(ix).alt_5prime_details.exons(:,end+1) = [idx_exon_l;idx_exon_r;idx_exon_alt];
      exons_is_alt([idx_exon_l;idx_exon_r;idx_exon_alt]) = exons_is_alt([idx_exon_l;idx_exon_r;idx_exon_alt]) + 1;
    end
  else
    genes(ix).alt_5prime = 0;
  end

  % alternative 3prime
  if ismember(ix,idx_alt_3prime)
    genes(ix).alt_3prime_details.exons = zeros(3,0);
    tix = find(idx_alt_3prime==ix);
    genes(ix).alt_3prime = length(tix);
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
      genes(ix).alt_3prime_details.exons(:,end+1) = [idx_exon_l;idx_exon_r;idx_exon_alt];
      exons_is_alt([idx_exon_l;idx_exon_r;idx_exon_alt]) = exons_is_alt([idx_exon_l;idx_exon_r;idx_exon_alt]) + 1;
    end
  else
    genes(ix).alt_3prime = 0;
  end

  % xor exons
  if ismember(ix,idx_xor_exons)
    tix = find(idx_xor_exons==ix);
    genes(ix).alt_xor_exon = length(tix);
    for ixt = tix
      exons_is_alt(exon_xor_exons(:,tix)) = exons_is_alt(exon_xor_exons(:,tix)) + 1;
    end
  else
    genes(ix).alt_xor_exon = 0;
  end

  % alternative transcriptions
  if ismember(ix,idx_alt_tstart)
    tix = find([idx_alt_tstart == ix]);
    % there is only one entry per gene.
    assert(all(size(tix) == [1,1]))
    genes(ix).alt_tstart = length(exon_alt_tstart{tix}.exons);
    genes(ix).alt_tstart_details.exons = exon_alt_tstart{tix}.exons;
  else
    genes(ix).alt_tstart = 0;
  end

  if ismember(ix,idx_alt_tend)
    tix = find([idx_alt_tend == ix]);
    % there is only one entry per gene.
    assert(all(size(tix) == [1,1]))
    genes(ix).alt_tend = length(exon_alt_tend{tix}.exons);
    genes(ix).alt_tend_details.exons = exon_alt_tend{tix}.exons;
  else
    genes(ix).alt_tend = 0;
  end

  genes(ix).exons_is_alt = exons_is_alt;
end


if 0,
  fid = fopen(sprintf('%s/alt_splice_type.tsv', genome_info.basedir),'w');
  fprintf(fid,'Clust\tExon Skip\tIntron Reten\tAlt 5prime\tAlt 3prime\n');
  for ix = idx_alt  
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\n',ix,...
	    genes(ix).alt_exon,genes(ix).alt_intron,genes(ix).alt_5prime,genes(ix).alt_3prime);
  end
  fclose(fid);
end
