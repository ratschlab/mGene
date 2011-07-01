function generate_listofNs(fn_in,fn_out,type)
%  generate_listofNs(fn_in,fn_out,type)
 
if nargin<3
  type = 'NSPACER';
end
  
tmp_fname = tempname;
unix(sprintf('grep %s %s > %s', type, fn_in, tmp_fname));

fd1 = fopen(tmp_fname,'r');
fd2 = fopen(fn_out,'w+');
source = '.';
score = '.';
strand = '.';
phase = '.';
attr_str = '.';
cnt = 0;
while ~feof(fd1),
  cnt = cnt+1;
  Line = fgetl(fd1);
  if isempty(Line),continue, end
  if ~ischar(Line), break, end ;
  elems = separate(Line);
  chr = elems{1};
  start = elems{2};
  stop = elems{3};
  fprintf(fd2,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
          chr,source,type,start,stop,score,strand, ...
          phase,attr_str);
end

fclose(fd1);
fclose(fd2);