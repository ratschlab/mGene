function [chr, strand, pos, score, content_name, score_name, chr_name_dic, content_name_dic, score_name_dic, success]=read_contpred_bin(contpred_fname)
% [chr, strand, pos, score, content_name, score_name, chr_name_dic, content_name_dic, score_name_dic, success]=read_contpred_bin(contpred_fname)

chr_name_dic = {} ;
score_name_dic = {} ;
content_name_dic = {} ;

chr=[] ; 
content_name = [] ;
score_name = [] ;
strand = '' ;
pos = [] ;
score = [] ;

if ~fexist(contpred_fname),
  success=false ;
  return ;
end ;

cnt = 0; last=0;
while ~last,
  try
    CPF=load(contpred_fname, sprintf('CPF_%i', cnt)) ;
    if ~isfield(CPF, sprintf('CPF_%i', cnt)),
      break
    end ;
  catch
    break ;
  end 
  if ~isempty(CPF)
    eval(sprintf('CPF=CPF.CPF_%i;', cnt)) ;
    cnt=cnt+1 ;
  %else
  %  CPF=load(contpred_fname) ;
  %  last = 1 ;
  end ;

  if isempty(CPF.pos)
    assert(isempty(CPF.score))
    continue;
  end

  if any(ismember(chr_name_dic, CPF.contig_name)),
    contig_idx = find(ismember(chr_name_dic, CPF.contig_name)) ;
  else
    chr_name_dic{end+1} = CPF.contig_name ;
    contig_idx = length(chr_name_dic) ;
  end ;

  if any(ismember(score_name_dic, CPF.score_name)),
    score_idx = find(ismember(score_name_dic, CPF.score_name)) ;
  else
    score_name_dic{end+1} = CPF.score_name ;
    score_idx = length(score_name_dic) ;
  end ;

  if any(ismember(content_name_dic, CPF.content_name)),
    content_idx = find(ismember(content_name_dic, CPF.content_name)) ;
  else
    content_name_dic{end+1} = CPF.content_name ;
    content_idx = length(content_name_dic) ;
  end ;

  chr(end+1:end+size(CPF.pos,1))=contig_idx ;
  score_name(end+1:end+size(CPF.pos,1))=score_idx ;
  content_name(end+1:end+size(CPF.pos,1))=content_idx ;
  pos(end+1:end+size(CPF.pos,1),:)=CPF.pos ;
  score(end+1:end+size(CPF.pos,1))=CPF.score ;
  strand(end+1:end+size(CPF.pos,1))=CPF.strand ;
end ;

success = true ;
