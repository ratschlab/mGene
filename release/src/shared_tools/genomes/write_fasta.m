function write_fasta(fd, name, sequence, end_line1, end_line2) ;

% function write_fasta(fd, name, sequence, end_line1, end_line2) ;
% 
% 
% INPUT - fd:                       file identifier 
%       - name:                     sequence name
%       - sequence
%       - end_line1 (default '\n'): following sequence name 
%       - end_line2 (default '\n') : following sequence
% 
% SEE ALSO  write_genes_cds2fasta 
  
if nargin<4, end_line1 = sprintf('\n') ; end ;
if nargin<5, end_line2 = sprintf('\n') ; end ;

fprintf(fd,'>%s%s', name, end_line1) ;
for i=1:80:length(sequence)
  fprintf(fd, '%s%s', sequence(i:min(length(sequence),i+79)), end_line2) ;
end ;
