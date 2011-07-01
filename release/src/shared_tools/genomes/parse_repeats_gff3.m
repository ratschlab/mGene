function genes = parse_repeats_gff3(gff_fname, save_fname,source)
% genes = parse_gff3(gff_fname, save_fname,source)
  
  
if nargin<3
  source = [];
end

%%% GFF file summary

[fd msg] = fopen(gff_fname, 'r');
disp(msg);

fprintf('\n------FILE SUMMARY-----\n');
[tmp,seqids] = unix(sprintf('sed -e "1d" %s | cut -f 1 |sort -u',gff_fname));
fprintf('gff file containes following seqids in column 1 :\n');
fprintf('%s\n',seqids)

[tmp,sources]=unix(sprintf('sed -e "1d" %s | cut -f 2 |sort -u',gff_fname));
fprintf('gff file containes following sources in column 2 :\n');
sources = separate(sources,sprintf('\n'));
sources(end)=[];
for t=1:length(sources)
  [tmp,num_sources]=unix(sprintf('cut -f 2 %s |grep %s | wc -l',gff_fname,sources{t}));
  fprintf('%s\t%i\n',sources{t},str2num(num_sources));
end
  
if ~isempty(source) 
  fprintf('filtering for source :\n');
  fprintf('%s\n',source);
  if all(isempty(strfind(sources,source)))
    fprintf('source not found in file\n')
    return
  end
  tmp_fname = tempname;
  unix(sprintf('grep %s %s > %s', source, gff_fname, tmp_fname));
else
  tmp_fname = gff_fname;
end

[tmp,types]=unix(sprintf('sed -e "1d" %s | cut -f 3 |sort -u',tmp_fname));
fprintf('\nfiltered gff file containes following types in column 3 :\n');
types = separate(types,sprintf('\n'));
types(end)=[];

for t=1:length(types)
  [tmp,num_types]=unix(sprintf('cut -f 3 %s | grep %s|wc -l',tmp_fname,types{t}));
  fprintf('%s\t%i\n',types{t},str2num(num_types));
end


if fexist(sprintf('%s.mat',save_fname))
  fprintf('file %s already exists\n',save_fname)
  fprintf('\nloaded repeat data from %s\n', save_fname);
  load(save_fname,'repeats');
  assert(length(repeats)==num_repeats) 
  return
end


fprintf('\n------starting parse\n');
for i=1:length(types)
  [tmp,num] = unix(sprintf('cut -f 3 %s | grep %s |wc -l',tmp_fname,types{i}));
  num_repeats(i) = str2num(num);
  fprintf('%s\t%i\n',types{i},num_repeats(i));
end
num_repeats = sum(num_repeats);

repeats(num_repeats).id = [];
repeats(num_repeats).name = [];
repeats(num_repeats).chr = [];
repeats(num_repeats).chr_num = [];
repeats(num_repeats).start = [];
repeats(num_repeats).stop = [];
repeats(num_repeats).score = [];
repeats(num_repeats).info = [];
[fd msg] = fopen(tmp_fname, 'r');
repeats_cnt = 0;
while ~feof(fd),
  parse_starter
  identifier = '';    
  for i=1:length(attributes),
    id=strfind(attributes{i},'ID=');
    if ~isempty(id),
      identifier = attributes{i}(id+3:end);
    end 
  end
  %   assert(~isempty(identifier)) ;
  
  repeats_cnt =  repeats_cnt +1;
  repeats(repeats_cnt).id = repeats_cnt;
  if ~isempty(identifier)
    repeats(repeats_cnt).name = identifier;
  end
  repeats(repeats_cnt).chr = elems{1};
  repeats(repeats_cnt).start = str2num(elems{4});
  repeats(repeats_cnt).stop = str2num(elems{5});
  repeats(repeats_cnt).score = str2num(elems{6});
  repeats(repeats_cnt).info = elems{9};
  fprintf('  parsed %i repeats\r', repeats_cnt);  
end % end first parse
% assert(repeats_cnt==num_repeats)
repeats=repeats(1:repeats_cnt);
fclose(fd);
fprintf('\n\n');
  
  
save(save_fname, 'repeats');


