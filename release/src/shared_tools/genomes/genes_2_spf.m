function [] = genes_2_spf(genes_fname,out_fname,signal_name)
% [] = genes_2_spf(genes_file,outfile,signal_name)
% generates a SPF file from a genes file

% test if filenames are set
if nargin < 3
  warning('to few parameters');
  return
end

% test if infile is existent
if isempty(genes_fname)
  warning('no infile specified');
  return;
elseif ~fexist(sprintf('%s',genes_fname))
  warning('file not existent or maybe wrong path');
  return;
end

% test if file 'out_fname' already exists
if isempty(out_fname)
  warning('no outfile specified');
  return;
elseif fexist(sprintf('%s',out_fname))
  warning(sprintf('file %s already exists. Please remove!\n',out_fname));
  return;
end

% test if signal_name is a valid signal name
if ~any(strcmp({'tss','tis','acc','don','cdsStop','cleave'}, signal_name))
  warning('wrong signal_name - possible choices include tss, tis, acc, don, cdsStop, cleave.');
end

% Description of SPF format taken from GALAXY server
% SPF format has exacty six required fields per line:

% 1. chrom - The name of the chromosome (e.g. chr1, chrY_random).
% 2. signalName - possible choices include tss, tis, acc, don, cdsStop, cleave.
% 3. scoreName - possible choices include label, output, Conf, Conf_Cum
% 4. chromPos - The position in the chromosome. (The first base in a chromosome is numbered 1.)
% 5. strand - Defines the strand - either '+' or '-'.
% 6. score - The score between -infinity and infinity. If scoreName is 'label', then the score should be -1/1 for negativ/postitiv.
cd 
load('-MAT',genes_fname,'genes');
% a field genes.config has to be added?
genes(1).config='/fml/ag-raetsch/share/databases/genomes/H_sapiens/hg18/genebuild/genome.config';
genes=load_sequence(genes);
fid=fopen(out_fname,'w');

for j=1:length(genes)
  if ~genes(j).is_valid
    continue;
  end
  % for all transcripts
  for i=1:length(genes(j).transcripts)
    % Build a line in output format
    pos_example=sprintf('%s\t%s\tlabel\t%u\t%s\t1.000\n',genes(j).chr, signal_name, genes(j).(signal_name)(i), genes(j).strand);
    % Write it to outfile
    fwrite(fid,pos_example);
    % search for negative examples in gene sequence
    pos=strfind(genes(j).seq ,'ATG');
    %length(pos)
    pos(find(pos == genes(j).tss(i)-genes(j).(signal_name)(i)+1))=[];
    %if strcmp(genes(j).strand,'+')
    %  pos(find(pos == genes(j).tis))=[];
    %else
    %  % need for translation of coordinates to '-' strand
    %  pos(find(pos == genes(j).tss-genes(j).tis))=[];
    %end
    %length(pos)
    for k=1:length(pos)
      neg_example=sprintf('%s\t%s\tlabel\t%u\t%s\t-1.000\n', genes(j).chr, signal_name, genes(j).tss(i)-pos(k)+1, genes(j).strand);
      fwrite(fid,neg_example);
    end
  end
end
fclose(fid);
fprintf('\ndone\n');
end