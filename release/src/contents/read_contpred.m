function [chr, strand, pos, score, content_name, score_name, chr_name_dic, content_name_dic, score_name_dic]=read_contpred(contpred_fname)
% [chr, strand, pos, score, content_name, score_name]=read_contpred(contpred_fname)

if fexist([contpred_fname '.gz']) && ~fexist(contpred_fname),
  ret = unix(sprintf('gzip -d %s',contpred_fname)) ;
  assert(ret==0) ;
end ;


chr_name_dic = {} ;
score_name_dic = {} ;
content_name_dic = {} ;

chr=[] ; 
content_name = [] ;
score_name = [] ;
strand = '' ;
pos = [] ;
score = [] ;

fd = fopen(contpred_fname, 'r') ;

rest_text = '' ; part = 0 ; part_offset = 0 ;
while ~feof(fd),
  part = part + 1 ;
  text=fread(fd, 100000000, 'uint8=>char') ;
  text=[rest_text text'] ;
  rest_text = '' ;
  if length(text)==0, 
    break ; 
  end ;
  last_idx = find(text==sprintf('\n'), 1, 'last') ;
  if last_idx~=length(text) && ~feof(fd),
    rest_text = text(last_idx+1:end) ;
    text(last_idx:end) = [] ;
  end ;
  text = [text sprintf('\n')] ;
  line_idx = find(text==sprintf('\n')) ;

  if length(line_idx)>2, 
    chr(length(chr)+1 : length(chr) + length(line_idx) -2 ) = nan ; 
    content_name(length(content_name)+1:length(content_name) + length(line_idx)-2) = nan ;
    score_name(length(score_name)+1 : length(score_name) + length(line_idx)-2) = nan ;
    strand(length(strand)+1 : length(strand) + length(line_idx)-2) = ' ' ;
    pos(length(pos)+1 : length(pos) + length(line_idx)-2,1:2) = nan ;
    score(length(score)+1 : length(score) + length(line_idx)-2) = nan ;
  end ;

  chr_idx = 0 ;
  content_name_idx = 0 ;
  score_name_idx = 0 ;
  progress = 0 ;
  for i=1:length(line_idx)-1,
    line=text(line_idx(i)+1:line_idx(i+1)-1) ;
    if isempty(line), 
      chr(i) = nan ;
      score(i) = nan ;
      score_name(i) = nan ;
      content_name(i) = nan ;
      pos(i,:) = nan ;
      strand(i) = ' ' ;
      continue; 
    end ;
    if isequal(line(1),'#')
      continue
    end
    elem_idx = find(line==sprintf('\t')) ;
    
    if (i-1)/length(line_idx)>=progress,
      fprintf('\rPart %i: Read line %i/%i (%2.0f%%)     ', part, (i-1), length(line_idx), progress*100);
      progress = progress+0.05 ;
    end ;
    if chr_idx==0 || ~isequal(line(1:elem_idx(1)-1), chr_name_dic{chr_idx}),
      chr_idx = strmatch(line(1:elem_idx(1)-1), chr_name_dic, 'exact') ;
      if isempty(chr_idx),
        chr_name_dic{end+1} = line(1:elem_idx(1)-1) ;
        chr_idx = length(chr_name_dic) ;
      end ;
    end ;
    
    if content_name_idx==0 || ~isequal(line(elem_idx(1)+1:elem_idx(2)-1), content_name_dic{content_name_idx}),
      content_name_idx = strmatch(line(elem_idx(1)+1:elem_idx(2)-1), content_name_dic, 'exact') ;
      if isempty(content_name_idx),
        content_name_dic{end+1} = line(elem_idx(1)+1:elem_idx(2)-1) ;
        content_name_idx = length(content_name_dic) ;
      end ;
    end ;
    
    if score_name_idx==0 || ~isequal(line(elem_idx(2)+1:elem_idx(3)-1), score_name_dic{score_name_idx}),
      score_name_idx = strmatch(line(elem_idx(2)+1:elem_idx(3)-1), score_name_dic, 'exact') ;
      if isempty(score_name_idx),
        score_name_dic{end+1} = line(elem_idx(2)+1:elem_idx(3)-1) ;
        score_name_idx = length(score_name_dic) ;
      end ;
    end ;
    
    chr(i+part_offset) = chr_idx ;
    content_name(i+part_offset) = content_name_idx ;
    score_name(i+part_offset) = score_name_idx ;
    strand(i+part_offset) = line(elem_idx(4)+1:elem_idx(5)-1) ;
    pos(i+part_offset,1) = str2double(line(elem_idx(3)+1:elem_idx(4)-1)) ;
    pos(i+part_offset,2) = str2double(line(elem_idx(4)+1:elem_idx(5)-1)) ;
    score(i+part_offset) = str2double(line(elem_idx(6)+1:end)) ;
  end

  part_offset = part_offset + length(line_idx)-1 ;
end ;
fclose(fd) ;

fprintf('\rRead line %i/%i (%2.0f%%)              \n', length(line_idx), length(line_idx), 100);
content_name(isnan(chr)) = [] ;
score_name(isnan(chr)) = [] ;
strand(isnan(chr)) = [] ;
pos(isnan(chr),:) = [] ;
score(isnan(chr)) = [] ;
chr(isnan(chr)) = [] ;

assert(all(chr~=0)) ;

