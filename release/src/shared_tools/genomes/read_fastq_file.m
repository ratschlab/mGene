function [strs,quals,names]=read_fastq_file(fname)
% [strs,quals,names]=read_fastq_file(fname)
% 
% reads all entries in a fastq file
  
fd=fopen(fname,'r') ;

strs={} ;
quals={} ;
names={} ;

while(1)
  [strs{end+1}, quals{end+1}, names{end+1},last] = read_fastq(fd) ;
  if isempty(strs{end})
    assert(last==1) ;
    strs(end) = [] ;
    quals(end) = [] ;
    names(end) = [] ;
  end ;
  if last==1, break ; end ;
  if feof(fd), break ; end ;
end ;

fclose(fd) ;
