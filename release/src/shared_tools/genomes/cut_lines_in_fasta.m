function cut_lines(fasta_fname,fasta_fname_1)
  

[fd msg] = fopen(fasta_fname, 'r');
while ~feof(fd),
  Line = fgetl(fd);
  if isequal(Line(1),'>')
    if exist('fd1')
      fclose(fd1);
    end
    fasta_fn = sprintf('%sCHROMOSOME_%s.dna',fasta_fname_1,Line(2:end));
    [fd1 msg] = fopen(fasta_fn, 'w');
  end
  len = length(Line);
  num=0;
  while num<len
    fprintf(fd1,'%s\n',Line(num+1:min(num+80,len)));
    num=num+80;
  end
end
fclose(fd);
fclose(fd1);

return  
  
[fd msg] = fopen(fasta_fname, 'r');
[fd1 msg] = fopen(fasta_fname_1, 'w');
while ~feof(fd),
  Line = fgetl(fd);
  len = length(Line);
  num=0;
  while num<len
    fprintf(fd1,'%s\n',Line(num+1:min(num+80,len)));
    num=num+80;
  end
end 
fclose(fd);
fclose(fd1);

  
return  