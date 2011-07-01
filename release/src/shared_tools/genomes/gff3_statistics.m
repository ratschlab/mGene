function [tmp_fname,STATS] = gff3_statistics(gff_fname,thesources, contig)
% [tmp_fname,STATS] = gff3_statistics(gff_fname,thesources, contig)
  

if nargin<2
  thesources = [];
end

                                                                                                                                                       
if nargin<3                                                                                     
  contig = [];                                                                                    
end 
                                                                                    

%%% SOFA mapping
SOFA_short_names = {'gene','mRNA','exon','CDS','intron','polyA_sequence','polyA_site','five_prime_UTR','three_prime_UTR'};
SOFA_IDs = {'SO:0000704','SO:0000234','SO:0000147','SO:0000316','SO:0000188','SO:0000610','SO:0000553','SO:0000204','SO:0000205'};

%%% GFF file summary

[fd msg] = fopen(gff_fname, 'r');
disp(msg);

tmp_fname = tempname;
if ~isempty(thesources) 
  fprintf('filtering for source :\n');
  unix(sprintf('echo > %s', tmp_fname)) ;
  for i=1:length(thesources),
    fprintf('%s\n',thesources{i});
    unix(sprintf('grep %s %s |grep -v "#" >> %s', thesources{i}, gff_fname, tmp_fname));
  end 
else
  unix(sprintf('grep -v "#" %s > %s', gff_fname, tmp_fname));
end

fprintf('\n------FILE SUMMARY-----\n');
[tmp,seqids] = unix(sprintf('sed -e "1d" %s | cut -f 1 |sort -u', tmp_fname));
if isempty(contig) 
   fprintf('gff file containes following %i seqids in column 1 :\n',length(seqids));
   if length(seqids)<100
      fprintf('%s\n',seqids)
   end
end

[tmp, sources]=unix(sprintf('sed -e "1d" %s | cut -f 2 |sort -u', tmp_fname));
if isempty(sources)
  error('specified sources not contained in file')
end
fprintf('gff file containes following sources in column 2 :\n');
sources = separate(sources,sprintf('\n'));
sources(end)=[];
for t=1:length(sources)
  [tmp,num_sources]=unix(sprintf('cut -f 2 %s |grep %s | wc -l',tmp_fname,sources{t}));
  fprintf('%s\t%i\n',sources{t},str2num(num_sources));
end
  
if ~isempty(contig) 
  cnamepos = findstr(upper(contig),upper(seqids));
  % contig = seqids(cnamepos:cnamepos+length(contig)-1); 
  fprintf('filtering for contig :\n');
  fprintf('^%s\n',contig);
  if all(isempty(strfind(upper(seqids),upper(contig))))
    fprintf('contig not found in file\n')
    return
  end
  tmp_fname1 = tempname;
  unix(sprintf('grep "^%s\t" %s > %s', contig, tmp_fname, tmp_fname1));
  [tmp,num_seqs] = unix(sprintf('sed -e "1d" %s | cut -f 1 |sort -u | wc -l',tmp_fname1));
  if (str2num(num_seqs)~=1)
    error('grep for contig doesnt work, repeat with by_contig=0')
  end
else
  tmp_fname1 = tmp_fname;
end
[tmp,types_]=unix(sprintf('sed -e "1d" %s | cut -f 3 |sort -u',tmp_fname1));
fprintf('\nfiltered gff file containes following types in column 3 :\n');
types = separate(types_,sprintf('\n')) ;
idx = find(cellfun('isempty',types));
types(idx)=[];


for t=1:length(types)
  [tmp,num_types]=unix(sprintf('cut -f 3 %s | grep ^%s\t |wc -l',tmp_fname1,types{t}));
  fprintf('%s\t%i\n',types{t},str2num(num_types));
  STATS.(types{t})=str2num(num_types);
end
% assert(all(ismember(lower(types),lower(SOFA_short_names)))|all(ismember(types,SOFA_IDs)))

tmp_fname= tmp_fname1;
