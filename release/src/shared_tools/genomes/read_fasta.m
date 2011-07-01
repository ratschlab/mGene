function [str,name,last]=read_fasta(fd) ;
% [str,name,last]=read_fasta(fd) ;
%
% reads the next entry in a fasta file

line='' ;
while isempty(line)
  line=fgetl(fd) ;
end 
if ~ischar(line),
  str='' ;
  name='' ;
  last=1 ;
  return ;
end ;

assert(line(1)=='>') ;
name=line(2:end) ;

str=char(ones(1,10000));
len = 1;
pos=ftell(fd) ;
line=fgetl(fd) ;
while ischar(line) &  ~isempty(line) & line(1)~='>',
  new_len = length(line);
  if len+new_len>length(str)
    str = [str char(ones(1,length(str)))];
  end
  str(len:len+new_len-1) = line;
  len = len+new_len;
  pos=ftell(fd) ;
  line=fgetl(fd) ;
end ;
str(len:end) = [];
fseek(fd,pos,-1) ;
if ~ischar(line),
  last=1 ;
else
  last=0 ;
end ;
