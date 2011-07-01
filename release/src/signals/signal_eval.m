function [ROCs, PRCs, signal_score, signal_score_label] = signal_eval(spf_label_fname, spf_label_dir, spf_pred_fname, spf_pred_dir, info_fname)
% [ROCs, PRCs] = galaxy_signal_eval(fn_label, fn_pred, info_fname)
%

disp('-----------------------------') ;
disp('Step 1: Reading label file...') ;
disp('-----------------------------') ;

spf_label_mat_fname = sprintf('%s/label_spf.mat', spf_label_dir) ;
if fexist(spf_label_mat_fname),
  disp('Reading binary SPF label files') ;
  [chr_label, strand_label, pos_label, score_label, score_names_label, signal_names_label, chr_label_dic, signal_names_label_dic, score_names_label_dic, success]=read_sigpred_bin(spf_label_mat_fname) ;
  assert(success) 
else
  disp('Reading SPF label files') ;
  [chr_label, strand_label, pos_label, score_label, score_names_label, signal_names_label, chr_label_dic, signal_names_label_dic, score_names_label_dic]=read_sigpred(spf_label_fname) ;
end ;

signal_names = unique(signal_names_label_dic) ;

idx_ = strmatch('label', score_names_label_dic, 'exact') ;
idx = find(score_names_label==idx_) ;
if length(idx)~=length(score_names_label),
  warning('dropping %i/%i non-label entries from label file', length(score_names_label)-length(idx), length(score_names_label)) ;
  chr_label = chr_label(idx) ;
  strand_label = strand_label(idx) ;
  pos_label = pos_label(idx) ;
  score_label = score_label(idx) ;
  signal_names_label = signal_names_label(idx) ;
  score_names_label = score_names_label(idx) ;
end 

idx = find(abs(abs(score_label)-1)<=1e-3) ;
if length(idx)~=length(score_label),
  warning('dropping %i/%i entries with label not +/-1 from label file', length(score_label)-length(idx), length(score_label)) ;
  chr_label = chr_label(idx) ;
  strand_label = strand_label(idx) ;
  pos_label = pos_label(idx) ;
  score_label = score_label(idx) ;
  signal_names_label = signal_names_label(idx) ;
  score_names_label = score_names_label(idx) ;
end 
score_label = sign(score_label) ;

fprintf('Done.\n\n') ;

disp('----------------------------------') ;
disp('Step 2: Reading prediction file...') ;
disp('----------------------------------') ;

spf_pred_mat_fname = sprintf('%s/Conf_cum_spf.mat', spf_pred_dir) ;
if fexist(spf_pred_mat_fname),
  disp('Reading binary SPF prediction files') ;
  [chr, strand, pos, score, signal_names, score_names, chr_dic, signal_names_dic, score_names_dic, success] = read_sigpred_bin(spf_pred_mat_fname) ;
  assert(success) ;
else
  disp('Reading SPF prediction files') ;
  [chr, strand, pos, score, signal_names, score_names, chr_dic, signal_names_dic, score_names_dic] = read_sigpred(spf_pred_fname) ;
end ;

eval_signal_names = intersect(signal_names_dic, signal_names_label_dic) ;

fprintf('Done.\n\n') ;

disp('---------------------------------') ;
disp('Step 3: Evaluating predictions...') ;
disp('---------------------------------') ;

% determine common chr numbering
[eval_chr, chr_map, chr_label_map] = intersect(lower(chr_dic), lower(chr_label_dic))  ;
chr = -chr ;
for i=1:length(chr_map),
  idx=find(-chr==chr_map(i)) ;
  chr(idx) = i ;
end ;
chr(chr<=0) = nan ;

chr_label = -chr_label ;
for i=1:length(chr_label_map),
  idx=find(-chr_label==chr_label_map(i)) ;
  chr_label(idx) = i ;
end ;
chr_label(chr_label<=0) = nan ;

fi = fopen(info_fname, 'w+') ;
ROCs=[] ;
PRCs=[] ;
for i=1:length(eval_signal_names),
  signal_name = eval_signal_names{i} 
  fprintf(fi, 'Evaluation for "%s" predictions ', signal_name) ;

  idx_label_ = strmatch(signal_name, signal_names_label_dic, 'exact') ;
  idx_label = find(signal_names_label==idx_label_ & ~isnan(chr_label)) ;
  idx_ = strmatch(signal_name, signal_names_dic, 'exact') ;
  idx = find(signal_names == idx_ & ~isnan(chr)) ;

  %id = [chr(idx); double(strand(idx)); pos(idx)]' ;
  %id_label = [chr_label(idx_label); double(strand_label(idx_label)); pos_label(idx_label)]' ;

  id_ = (chr(idx)*100000000+pos(idx)).*(1-2*(strand(idx)=='-')) ;
  id_label_ = (chr_label(idx_label)*100000000+pos_label(idx_label)).*(1-2*(strand_label(idx_label)=='-')) ;

  [tmp, idx_id, idx_id_label] = intersect(id_, id_label_);%, 'rows') ;
  fprintf(fi, 'on %i positions\n', length(idx_id)) ;

  signal_score = score(idx(idx_id)) ;
  signal_score_label = score_label(idx_label(idx_id_label)) ;

  pos_ = pos(idx(idx_id)) ;
  pos_label_ = pos_label(idx_label(idx_id_label)) ;
  assert(isequal(pos_, pos_label_)) ;

  ROC=calcrocscore(signal_score, signal_score_label) 
  PRC=calcrfcscore(signal_score, signal_score_label) 

  fprintf(fi, ' * Area under ROC curve: %1.3f\n', ROC) ;
  fprintf(fi, ' * Area under PRC curve: %1.3f\n\n', PRC) ;

  ROCs = [ROCs ROC] ;
  PRCs = [PRCs PRC] ;
end ;

fclose(fi) ;

fprintf('Done.\n\n') ;
