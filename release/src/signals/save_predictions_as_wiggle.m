function fn_out = save_predictions_as_wiggle(out_fn_prefix, pos, score, score_name, signal_name, contig_name, strand, do_gzip)
% fn_out = save_predictions_as_wiggle(out_fn_prefix, pos, score, score_name, signal_name, contig_name, strand, do_gzip)

%assert(issorted(pos))
assert(isequal(pos, sort(pos)))
assert(length(pos)==length(score))

if ~exist('do_gzip', 'var'),
  do_gzip = 1 ;
end ;

fn_out=[out_fn_prefix '_' score_name '.wiggle'];

if fexist(fn_out)|fexist([fn_out '.gz'])
  fprintf(1, 'file exists, it will NOT be overwritten\n')
%  keyboard
  return
end

fprintf(1,'save as wiggle: %s\n', fn_out);
[fd msg] = fopen(fn_out, 'w+');
disp(msg)

fprintf(fd,'browser position chr%s:1-50,000\n', contig_name);
fprintf(fd,'browser hide all\n');
fprintf(fd,'browser pack Genes\n');
fprintf(fd,'browser full altGraph\n');
fprintf(fd,'browser viewLimits=%1.2f:%1.2f\n', min(score), max(score));


fprintf(fd,['track type=wiggle_0 name="%s_%s%s" description="%s chr%s %s" ' ...
            '\n'],signal_name, score_name, strand, score_name, ...
        contig_name, strand);
fprintf(fd,'variableStep chrom=chr%s\n',contig_name);

for j=1:length(pos)
  fprintf(fd,'%i\t%.3f\n', pos(j), score(j));
end

fclose(fd);

if do_gzip,
  ret = unix(sprintf('gzip %s',fn_out));
  if ret~=0
    warning('gzip failed')
  end
end ;
