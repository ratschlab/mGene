function save_predictions_as_contpred(out_fn_prefix, pos, score, score_name, content_name, contig_name, strand, do_gzip)
% save_predictions_as_wiggle(out_fn_prefix, pos, score, score_name, signal_name, contig_name, strand, do_gzip)

%assert(issorted(pos))
  
assert(isempty(pos) || isequal(pos(:,1), sort(pos(:,1))))
assert(size(pos,1)==length(score))

if ~exist('do_gzip', 'var'),
  do_gzip = 1 ;
end ;

fn_out=[out_fn_prefix '_' score_name '.cpf'];

if fexist(fn_out) | fexist([fn_out '.gz'])
  fprintf(1, 'file exists, it will NOT be overwritten\n')
  return
end
fprintf(1,'save in content prediction format (cpf): %s\n', fn_out);
[fd msg] = fopen(fn_out, 'w+');

fprintf(fd, '##Seqid\tType\tSource\tStart\tEnd\tStrand\tScore\n') ;

for j=1:min(10000,size(pos, 1))
  fprintf(fd, '%s\t%s\t%s\t%i\t%i\t%c\t%.3f\n', contig_name, content_name, score_name, pos(j,1),pos(j,2), strand, score(j)) ;
end
if length(pos)>10000,
  fprintf(fd, '##Attention: this file has been truncated.\n') ;
end ;
fclose(fd);

if do_gzip,
  ret = unix(sprintf('gzip %s',fn_out)) ;
  if ret~=0
    warning('gzip failed') ;
  end
end 

