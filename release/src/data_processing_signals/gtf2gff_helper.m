function gtf2gff(gtf_fname, gff3_fname) ;
% gtf2gff(gtf_fname, gff3_fname) ;

tmpname = tempname ;
[res] = unix(sprintf('cut -f 9 %s | sort -u> %s', gtf_fname, tmpname));
assert(res==0) ;

gene_ids = {} ;
gene_info = struct ;
gene_info.chr='' ;
gene_info.start=inf ;
gene_info.stop=0 ;
gene_info.strand='' ;
gene_info.info='' ;
gene_info=gene_info([]) ;

transcript_ids = {} ;
fd=fopen(gtf_fname, 'r') ;
fd2=fopen(gff3_fname, 'w') ;
while (~feof(fd)),
  line=fgetl(fd) ;
  if ~ischar(line),
    break ;
  end ;
  if line(1)=='#', continue; end ;
  items1=separate(line) ;
  chr=items1{1} ;
  info=items1{2} ;
  start=str2num(items1{4}) ;
  stop=str2num(items1{5}) ;
  strand=items1{7} ;
  
  items=separate(items1{9}, ';') ;
  gene_id='' ;
  for i=1:length(items),
    items{i}=deblank(items{i}) ;
    gene_id_pos = strfind(items{i}, 'gene_id') ;
    if ~isempty(gene_id_pos),
      str=items{i}(gene_id_pos+7:end) ;
      str(str=='"') = [] ;
      while(str(1)==' ')
        str=str(2:end) ;
      end ;
      gene_id=deblank(str) ;
      if length(gene_ids)>1 && isequal(gene_ids{end}, gene_id),
        idx=length(gene_ids) ;
      else
        idx=strmatch(gene_id, gene_ids, 'exact') ;
      end ;
      if isempty(idx)
        gene_ids{end+1}=str ;
        gene_info(length(gene_ids)).chr=chr ;
        gene_info(length(gene_ids)).start=start ;
        gene_info(length(gene_ids)).stop=stop ;
        gene_info(length(gene_ids)).strand=strand ;
        gene_info(length(gene_ids)).info=info ;
        transcript_ids{end+1}={} ;
        idx=length(gene_ids) ;
        if mod(idx,1000)==0, 
	  fprintf('reading gene %i\n', idx) ;
        end ;
      else
        gene_info(idx).start=min(start, gene_info(idx).start) ;
        gene_info(idx).stop=max(stop, gene_info(idx).stop) ;
      end ;
    end ;
  end ;
  if isempty(gene_id),
    fprintf('failed parsing line:\n%s\n',line) ;
    %save ~/tmp/gtf2gff_bug1.mat 
    assert(false);
  end ;
  transcript_id='' ;
  for i=1:length(items),
    items{i}=deblank(items{i}) ;
    transcript_id_pos = strfind(items{i}, 'transcript_id') ;
    if ~isempty(transcript_id_pos),
      str=items{i}(transcript_id_pos+13:end) ;
      str(str=='"') = [] ;
      while(str(1)==' ')
        str=str(2:end) ;
      end ;
      transcript_id = deblank(str) ;
      assert(length(idx)==1)
      if isempty(strmatch(transcript_id, transcript_ids{idx}))
        transcript_ids{idx}{end+1}=transcript_id ;
      end ;
    end ;
  end ;
  if isempty(transcript_id),
    printf('failed parsing line:\n%s\n',line) ;
    %save ~/tmp/gtf2gff_bug2.mat 
    assert(false);
  end ;
  
  fprintf(fd2, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tParent=transcript:%s\n', items1{1}, items1{2}, items1{3}, items1{4}, items1{5}, items1{6}, items1{7}, items1{8}, transcript_id);
end ;
fclose(fd) ;

for i=1:length(gene_ids)
  fprintf(fd2, '%s\t%s\tgene\t%i\t%i\t.\t%c\t.\tID=gene:%s\n', gene_info(i).chr, gene_info(i).info, gene_info(i).start, gene_info(i).stop, gene_info(i).strand, gene_ids{i}) ;
  for j=1:length(transcript_ids{i}),
    fprintf(fd2, '%s\t%s\tmRNA\t%i\t%i\t.\t%c\t.\tID=transcript:%s;Parent=gene:%s\n', gene_info(i).chr, gene_info(i).info, gene_info(i).start, gene_info(i).stop, gene_info(i).strand, transcript_ids{i}{j}, gene_ids{i}) ;
  end ;
end
fclose(fd2) ;

% cleanup
res=unix(sprintf('rm %s', tmpname)) ;
assert(res==0) ;

return

function x = fexist(filename,gz)
% FEXIST Check if file exists
%        FEXIST('filename') returns:
%        0 if filename does not exist
%        1 if filename exists in the given path

if nargin<2,
  gz=0;
end ;

fd = fopen(filename);
if fd ~= -1
  x = 1;
  fclose(fd);
elseif gz,
  fd = fopen([filename '.gz']);
  if fd~=-1,
    x=1 ;
    fclose(fd) ;
  else
    x=0 ;
  end ;
else
  x = 0;
end
x=logical(x) ;
return

