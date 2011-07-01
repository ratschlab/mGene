function [dname,fname] = find_internal_files(dname, fname) ;
% [dname,fname] = find_internal_files(dname, fname) ;

items = separate(dname, ':') ;

if (length(items)~=3 && length(items)~=2) || (length(dname)>1 && dname(1)=='/'),
  % leave things untouched
  return ;
end ;

switch items{1}
 case 'mgene'
 info_file = 'tools/mGene/internal/mgene_index.txt' ;
 case 'gene'
 info_file = 'tools/mGene/internal/gene_index.txt' ;
 case 'signal'
 info_file = 'tools/mGene/internal/signal_index.txt' ;
 case 'content'
 info_file = 'tools/mGene/internal/content_index.txt' ;
 case 'genome'
 info_file = 'tools/mGene/internal/genome_index.txt' ;
otherwise
 error('unknown token (%s)', items{1}) ;
end ;

fprintf('searching for token %s in index file: %s\n', dname, info_file) ;

fd = fopen(info_file, 'r') ;
if fd<1, error('could not open file %s', info_file); end
found = 0 ;
while (~feof(fd)),
  line = fgetl(fd);
  if ~ischar(line), break ; end ;
  line_items=separate(line) ;
  token_items=separate(line_items{1}, ':') ;
  assert(length(token_items)+1==length(items)) ;

  if length(token_items)==1 && isequal(items{2}, token_items{1}),
    dname = line_items{2} ;
    if length(line_items)>=3,
      fname = line_items{3} ;
    end ;
    found=1;
    break ;
  elseif length(token_items)==2 && isequal(items{2}, token_items{1}) && isequal(items{3}, token_items{2})
    dname = line_items{2} ;
    if length(line_items)>=3,
      fname = line_items{3} ;
    end ;
    found=1;
    break ;
  end ;
end ;

if ~isempty(dname) && dname(1)~='/',
  dname = [pwd '/' dname] ;
end ;
if ~isempty(fname) && fname(1)~='/',
  fname = [pwd '/' fname] ;
end ;

%[engine, environment]=determine_engine() ;
%if isequal(environment, 'internal'),
%  dname = strrep(dname, '/home/galaxy/internal_data', '/fml/ag-raetsch/share/software/galaxy-110509/internal_data') ;
%  fname = strrep(fname, '/home/galaxy/internal_data', '/fml/ag-raetsch/share/software/galaxy-110509/internal_data') ;
%end ;

if ~found,
  fprintf(2, 'internal data source not found (%s)', dname) ;
  exit(-1) ;
else
  fprintf('Located internal resources:\n\tresource1=%s\n\tresource2=%s\n', dname, fname) ;
end ;


