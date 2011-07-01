function [origisoforms,isoforms,isoORFlength,isostartpos,isostoppos] = find_all_orfs(isoforms,chr,strand,vertices)
% function [origisoforms,isoforms,isoORFlength,isostartpos,isostoppos] = find_all_orfs(isoforms,chr,strand,vertices)
%
% For each isoform, find the longest open reading frame.
% Sort by length of ORF, longest first.
% origisoforms contains all the exons, isoforms contains only the exons with coding sequence.

origisoforms = isoforms;
isostartpos = zeros(length(isoforms),1);
isostoppos = zeros(length(isoforms),1);
isoORFlength = zeros(length(isoforms),1);

for ixf = 1:length(isoforms)
  % find the longest open reading frames
  str = load_mrna(chr,strand,vertices(1,isoforms{ixf})',vertices(2,isoforms{ixf})');
  [strstart,strstop,isoORFlength(ixf)] = get_longest_orf(str,strand);

  % convert string start and stop to genomic locations
  [isostartpos(ixf),isostoppos(ixf),isoforms{ixf}] = strloc2genloc(isoforms{ixf}, vertices,...
                                                  strstart,strstop);
end

[dummy, iso_idx] = sort(isoORFlength,'descend');
isoORFlength = isoORFlength(iso_idx);
isostartpos = isostartpos(iso_idx);
isostoppos = isostoppos(iso_idx);
isoforms = isoforms(iso_idx);
origisoforms = origisoforms(iso_idx);
