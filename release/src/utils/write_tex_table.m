function write_tex_table(A, fd, col_labels, row_labels)

if nargin==2
  col_labels = repmat({'xxx'}, 1, size(A, 2));
  row_labels = repmat({'yyy'}, 1, size(A, 1));
end

for j = 1:length(col_labels)
	col_labels{j} = strrep(col_labels{j}, '_', '\_');
end
for j = 1:length(row_labels)
	row_labels{j} = strrep(row_labels{j}, '_', '\_');
end

rows = size(A, 1);
cols = size(A, 2);


include_header = 1;

if include_header
	fprintf(fd, '\\documentclass[a4paper, 12pt]{article}\n');
	fprintf(fd, '\\usepackage{times, mathptm}\n');
	fprintf(fd, '\\usepackage[latin1]{inputenc}\n');
	fprintf(fd, '\n');
	fprintf(fd, '\\begin{document}\n');
end



fprintf(fd, '\\begin{table}[hbt!]\n');
fprintf(fd, '  \\begin{center}\n');
fprintf(fd, '    \\caption{\\label{tab:xy}}\n');
fprintf(fd, '    \\par\n');
%fprintf(fd, '    \\scalebox{0.99}{\\mbox{\n');
fprintf(fd, '    \\begin{tabular}{|');
fprintf(fd, 'l|');
for c=2:cols+1, 
  fprintf(fd, 'r|');
end
fprintf(fd, '}\\hline\n');

fprintf(fd, '    &  \\multicolumn{%i}{c|}{\\bf title}\\\\\n',cols);
for c=1:cols,
  fprintf(fd, '    & {\\bf %s}', col_labels{c});
end
fprintf(fd, '\\\\\\hline\n');
for r=1:rows,
  fprintf(fd, '      %s ', row_labels{r}); 
  for c=1:cols
    fprintf(fd, '& %.1f ', A(r,c)*100);
    %fprintf(fd, '& %i ', A(r,c));
  end
  fprintf(fd, '\\\\\n');
end
fprintf(fd, '    \\hline\n');
fprintf(fd, '    \\end{tabular}\n');
%fprintf(fd, '    }}\n');
fprintf(fd, '  \\end{center}\n');
fprintf(fd, '\\end{table}\n');

if include_header
	fprintf(fd, '\n');
	fprintf(fd, '\\end{document}\n');
end

