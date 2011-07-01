function [info_names,sig,alt_names] = initialize_cand_search(Signal,info,genes,strand,signal_name,debug)
% [Signal,info,info_names,sig,alt_names] = initialize_cand_search(PAR,genes,strand,signal_name,debug)

rand('seed',1234);

% info = PAR.info_genes.(signal_name);
% Signal = PAR.Signals.(PAR.Signal_name);


info_names = fieldnames(info);
use = [];
fprintf(1,'getting conf for ');
for i=1:length(info_names)
  if ismember([info_names{i} '_conf'],Signal.Conf_names)
    use = [use i];
    if isequal(signal_name,'transcript')&~isequal(Signal.name,'transacc')
      assert(ismember([info_names{i} '_conf_skip'],Signal.Conf_names))
    end
    fprintf(1,'%s, ',info_names{i})
  end
end
fprintf('\n')
info_names = info_names(use);

% if ismember('cDNA',info_names)
%   use = strmatch('cDNA',info_names,'exact');
%   info_names{use,2} = 'fulllength';
%   fprintf(1,'\nfulllength cDNAs added to cDNAs\n')
% end

fprintf(1,'\n');

if ~isempty(genes)
  strands = unique([genes.strand]);
  assert(strand==strands(1))
  assert(length(strands)==1);
  assert(length([genes.chr_num])==length(genes));%all chr_num fields not empty
  assert(length(unique([genes.chr_num]))==1);
  if ~(length(unique([genes.chr_num]))==1)
	unique([genes.chr_num])
	error('genes not from a single chromosome')
  end
end


sig='';
switch Signal.name
 case {'acc','transacc','alt_exon_acc'}
  sig = 'acc';
 case {'don','alt_exon_don'}
  sig = 'don';
 otherwise
  sig = 'other';
end

if isequal(sig,'other') || isequal(Signal.name,'transacc') || ...
      isequal(Signal.name,'acc') || isequal(Signal.name,'don') 
  alt_names ={};
  return
end

alt_names = {'alt_exon', 'alt_intron', 'alt_5prime', 'alt_3prime'};
use = []; 
fprintf(1,'getting label 2 for ');
for i=1:length(alt_names)
  if ismember(alt_names{i},Signal.Conf_names)
    fprintf(1,'%s, ',alt_names{i})
  end
end
