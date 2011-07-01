
function hits = find_patterns( seq, consensi );

hits = [];
for( k = 1:length(consensi) )
  hits = [ hits , strfind(seq,consensi{k}) ];
end;
hits = sort( hits' );

