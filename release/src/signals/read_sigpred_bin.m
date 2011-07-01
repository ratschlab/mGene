function [chr, strand, pos, score, signal_name, score_name, chr_name_dic, signal_name_dic, score_name_dic, success]=read_sigpred_bin(sigpred_fname)
% [chr, strand, pos, score, signal_name, score_name, chr_name_dic, signal_name_dic, score_name_dic, success]=read_sigpred_bin(sigpred_fname)

chr_name_dic = {} ;
score_name_dic = {} ;
signal_name_dic = {} ;

chr=[] ; 
signal_name = [] ;
score_name = [] ;
strand = '' ;
pos = [] ;
score = [] ;

if ~fexist(sigpred_fname),
  success=false ;
  return ;
end ;

cnt = 0; last=0;
while ~last,
  try
    SPF=load(sigpred_fname, sprintf('SPF_%i', cnt)) ;
    if ~isfield(SPF, sprintf('SPF_%i', cnt)),
      break
    end ;
  catch
    break ;
  end 
  if ~isempty(SPF)
    eval(sprintf('SPF=SPF.SPF_%i;', cnt)) ;
    cnt=cnt+1 ;
  %else
  %  SPF=load(sigpred_fname) ;
  %  last = 1 ;
  end ;

  if any(ismember(chr_name_dic, SPF.contig_name)),
    contig_idx = find(ismember(chr_name_dic, SPF.contig_name)) ;
  else
    chr_name_dic{end+1} = SPF.contig_name ;
    contig_idx = length(chr_name_dic) ;
  end ;

  if any(ismember(score_name_dic, SPF.score_name)),
    score_idx = find(ismember(score_name_dic, SPF.score_name)) ;
  else
    score_name_dic{end+1} = SPF.score_name ;
    score_idx = length(score_name_dic) ;
  end ;

  if any(ismember(signal_name_dic, SPF.signal_name)),
    signal_idx = find(ismember(signal_name_dic, SPF.signal_name)) ;
  else
    signal_name_dic{end+1} = SPF.signal_name ;
    signal_idx = length(signal_name_dic) ;
  end ;

  chr(end+1:end+length(SPF.pos))=contig_idx ;
  score_name(end+1:end+length(SPF.pos))=score_idx ;
  signal_name(end+1:end+length(SPF.pos))=signal_idx ;
  pos(end+1:end+length(SPF.pos))=SPF.pos ;
  score(end+1:end+length(SPF.pos))=SPF.score ;
  strand(end+1:end+length(SPF.pos))=SPF.strand ;
end ;

success = true ;
