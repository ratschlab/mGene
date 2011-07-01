function save_predictions_as_sigpred_bin(out_fn_prefix, pos, score, score_name, signal_name, contig_name, strand)
% save_predictions_as_sigpred_bin(out_fn_prefix, pos, score, score_name, signal_name, contig_name, strand)

%assert(issorted(pos))
assert(isequal(pos, sort(pos)))
assert(length(pos)==length(score))

fn_out=[out_fn_prefix '_' score_name '_spf.mat'];

if fexist(fn_out) 
  fprintf(1, 'file exists, it will NOT be overwritten\n')
  return
end
fprintf(1,'save in binary signal prediction format (spf): %s\n', fn_out);
save(fn_out, '-V7', 'contig_name', 'signal_name', 'score_name', 'pos', 'strand', 'score') ;
