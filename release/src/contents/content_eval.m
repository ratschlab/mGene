function [ROCs, PRCs, score, LT] = content_eval(cpf_label_fname, cpf_label_dir, cpf_pred_fname, cpf_pred_dir, info_fname)
% [ROCs, PRCs] = content_eval(cpf_label_fname, cpf_label_dir, cpf_pred_fname, cpf_pred_dir, info_fname)
%

[ret, timedate] = unix('date') ;
fprintf('started content_eval %s', timedate) ;

fprintf('Results are written to: %s\n', info_fname) ;

disp('-----------------------------') ;
disp('Step 1: Reading label file...') ;
disp('-----------------------------') ;

cpf_label_mat_fname = sprintf('%s/label_cpf.mat', cpf_label_dir) ;
if fexist(cpf_label_mat_fname),
  disp('Reading binary CPF label files') ;
  [L.chr, L.strand, L.pos, L.score, L.score_names, L.content_names, ...
   L.chr_dic, L.content_names_dic, L.score_names_dic, success] = read_contpred_bin(cpf_label_mat_fname) ;
  assert(success) 
else
  disp('Reading CPF label files') ;
  [L.chr, L.strand, L.pos, L.score, L.score_names, L.content_names, ...
   L.chr_dic, L.content_names_dic, L.score_names_dic] = read_contpred(cpf_label_fname) ;
end ;

content_names = unique(L.content_names_dic) ;

idx_ = strmatch('label', L.score_names_dic, 'exact') ;
idx = find(L.score_names==idx_) ;
if length(idx)~=length(L.score_names),
  warning('dropping %i/%i non-label entries from label file', length(L.score_names)-length(idx), length(L.score_names)) ;
  L.chr = L.chr(idx) ;
  L.strand = L.strand(idx) ;
  L.pos = L.pos(idx,:) ;
  L.score = L.score(idx) ;
  L.content_names = L.content_names(idx) ;
  L.score_names = L.score_names(idx) ;
end 

idx = find(abs(abs(L.score)-1)<=1e-3) ;
if length(idx)~=length(L.score),
  warning('dropping %i/%i entries with label not +/-1 from label file', length(L.score)-length(idx), length(L.score)) ;
  L.chr = L.chr(idx) ;
  L.strand = L.strand(idx) ;
  L.pos = L.pos(idx,:) ;
  L.score = L.score(idx) ;
  L.content_names = L.content_names(idx) ;
  L.score_names = L.score_names(idx) ;
end 
L.score = sign(L.score) ;

fprintf('Done.\n\n') ;



disp('---------------------------------') ;
disp('Step 2: Evaluating predictions...') ;
disp('---------------------------------') ;

CHR =  unique(L.chr);
fi = fopen(info_fname, 'w+') ;
ROCs=[] ;
PRCs=[] ;

for i=1:length(L.content_names_dic),
  fprintf(fi, 'Evaluation for "%s" predictions ', L.content_names_dic{i}) ;
  spf_mat_fname = sprintf('%s/output_spf.mat', cpf_pred_dir) ;
  load(spf_mat_fname, 'SPF_info') ;

  score=[];
  LT=[];
  for c= CHR
    for s='+-'
      index = find(ismember({SPF_info(:).contig_name}, L.chr_dic{c}) & ([SPF_info(:).strand]==s)) 
      assert(length(index)==1) ;

      P=load(spf_mat_fname, sprintf('SPF_%i', index-1)) ;
      P=P.(sprintf('SPF_%i', index-1)) ;

      if ~isequal(P.signal_name,L.content_names_dic{i})||~isequal(P.contig_name,L.chr_dic{c})
        continue
      end

      idx = find(L.chr == c & L.strand == s);

      [tmp,idx_l1,idx_p1] = intersect(L.pos(idx,1),P.pos);
      [tmp,idx_l2,idx_p2] = intersect(L.pos(idx,2),P.pos);
      [tmp,ii1,ii2] = intersect(idx_l1,idx_l2);
      LT = [LT L.score(idx(idx_l1(ii1)))];
      if s=='+'  
        score = [score P.score(idx_p2(ii2)) - P.score(idx_p1(ii1))] ;
        calcrocscore(P.score(idx_p2(ii2)) - P.score(idx_p1(ii1)),L.score(idx(idx_l1(ii1)))) 
      else
        score = [score P.score(idx_p1(ii1)) - P.score(idx_p2(ii2))] ;
        calcrocscore(P.score(idx_p1(ii1)) - P.score(idx_p2(ii2)),L.score(idx(idx_l1(ii1)))) 
      end
    end
  end
  fprintf(fi, 'on %i positions (%i positives) \n', length(LT),sum(LT==1)) ;
  ROC = calcrocscore(score, LT) ; 
  PRC = calcrfcscore(score, LT) ;

  fprintf(fi, ' * Area under ROC curve: %1.3f\n', ROC) ;
  fprintf(fi, ' * Area under PRC curve: %1.3f\n\n', PRC) ;

  ROCs = [ROCs ROC] ;
  PRCs = [PRCs PRC] ;
end

fclose(fi) ;

fprintf('Done.\n\n') ;

[ret, timedate] = unix('date') ;
fprintf('finished content_eval %s', timedate) ;
