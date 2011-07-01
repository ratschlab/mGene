function [ee,ei,e5,e3] = find_unique_event(genes)
% function [ee,ei,e5,e3] = find_unique_event(genes)
%
% Find genes with only one type of alternative splicing event.

if 0,
  info = '~/altsplicedata/elegans/genome.config';
  origpath = '/fml/ag-raetsch/home/ong/splicing';

  addpath(sprintf('%s/utils', origpath)) ;
  addpath(sprintf('%s/sensors/altsplice', origpath)) ;
  addpath(sprintf('%s/splicegraphs', origpath)) ;
  addpath(sprintf('%s/splicegraphs/detect_altsplice', origpath)) ;
  addpath /fml/ag-raetsch/share/software/matlab_tools/utils
  
  genome_info = init_genome(info) ;
  fprintf('directory: %s\n', genome_info.basedir) ;
  load(sprintf('%s/confirmed_sequences.mat', genome_info.basedir), 'genes') ;
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

% find other alternatives which we don't consider now
[idx_alt_tstart, exon_alt_tstart] = detect_alttstart(genes) ;
[idx_alt_tend, exon_alt_tend] = detect_alttend(genes) ;
idx_alt_trans = union([idx_alt_tstart{:}],[idx_alt_tend{:}]);
[idx_xor_exons, exon_xor_exons] = detect_xorexons(genes, idx_alt) ;
[idx_multiple_skips, exon_multiple_skips] = detect_multipleskips(genes, idx_alt) ;
[idx_exon_intron5, exon_exon_intron5, idx_exon_intron3, exon_exon_intron3] = detect_exonintron(genes, idx_alt) ; 
[idx_alt_intron, introns_alt_intron] = detect_altintrons(genes, [1:length(genes)]) ;
idx_alt_other = union(union(union(union(idx_xor_exons,idx_multiple_skips),...
				  idx_exon_intron5),idx_exon_intron3),idx_alt_intron);

% find genes with only one alternative type
tee = setdiff(idx_exon_skips,union(union(idx_alt_5prime,idx_alt_3prime),idx_intron_reten));
tei = setdiff(idx_intron_reten,union(union(idx_alt_5prime,idx_alt_3prime),idx_exon_skips));
te3 = setdiff(idx_alt_3prime,union(union(idx_alt_5prime,idx_exon_skips),idx_intron_reten));
te5 = setdiff(idx_alt_5prime,union(union(idx_exon_skips,idx_alt_3prime),idx_intron_reten));

% remove the alternative transcription starts and stops
ee = setdiff(tee,idx_alt_trans);
ei = setdiff(tei,idx_alt_trans);
e3 = setdiff(te3,idx_alt_trans);
e5 = setdiff(te5,idx_alt_trans);

% remove the alternatives that we don't consider
ee = setdiff(ee,idx_alt_other);
ei = setdiff(ei,idx_alt_other);
e3 = setdiff(e3,idx_alt_other);
e5 = setdiff(e5,idx_alt_other);

