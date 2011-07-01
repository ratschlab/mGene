function save_predictions_as_contpred_bin(out_fn_prefix, pos, score, score_name, content_name, contig_name, strand)
% save_predictions_as_sigpred_bin(out_fn_prefix, pos, score, score_name, signal_name, contig_name, strand)

%assert(issorted(pos))

assert(isempty(pos) || isequal(pos(:,1), sort(pos(:,1))))
assert(size(pos,1)==length(score))

fn_out=[out_fn_prefix '_' score_name '_cpf.mat'];

if fexist(fn_out) 
  fprintf(1, 'file exists, it will NOT be overwritten\n')
  return
end
fprintf(1,'save in binary content prediction format (cpf): %s\n', fn_out);
save(fn_out, 'contig_name', 'content_name', 'score_name', 'pos', 'strand', 'score') ;
