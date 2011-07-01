function [str, qual, name, last]=read_fastq(fd) ;
% [str, qual, name, last]=read_fastq(fd) ;
%
% reads the next entry in a fastq file

line='' ;
while isempty(line)
  line=fgetl(fd) ;
end 
if ~ischar(line),
  str='' ;
  name='' ;
  qual='' ;
  last=1 ;
  return ;
end ;

assert(line(1)=='@') ;
name=line(2:end) ;

str=char(ones(1,10000));
len = 1;
pos=ftell(fd) ;
line=fgetl(fd) ;
num_lines = 0 ;
while ischar(line) && ~isempty(line) && line(1)~='+',
  new_len = length(line);
  if len+new_len > length(str)
    str = [str char(ones(1, length(str)))];
  end
  str(len:len+new_len-1) = line;
  len = len+new_len;
  line=fgetl(fd) ;
  num_lines = num_lines+1 ; ;
end ;
str(len:end) = [];

assert(line(1)=='+') ;
name2=line(2:end) ;
assert(isequal(name, name2)) ;


qual=char(ones(1,10000));
len = 1;
pos=ftell(fd) ;
line=fgetl(fd) ;
num_lines2 = 0 ;
while ischar(line) &&  ~isempty(line) && num_lines~=num_lines2%line(1)~='+' && line(1)~='@',
  new_len = length(line);
  if len+new_len>length(qual)
    qual = [qual char(ones(1,length(qual)))];
  end
  qual(len:len+new_len-1) = line;
  len = len+new_len;
  pos=ftell(fd) ;
  line=fgetl(fd) ;
  num_lines2=num_lines2+1 ;
end ;
qual(len:end) = [];

%assert(length(qual)==length(str)) ;

fseek(fd,pos,-1) ;
if ~ischar(line),
  last=1 ;
else
  last=0 ;
end ;

