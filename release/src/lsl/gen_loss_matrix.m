function truth = gen_loss_matrix(truth, all_pos, gene_confirmed, model)
% truth = gen_loss_matrix(truth, all_pos, gene_confirmed, model)
% gene_confirmed: -1: predicted
%                  0: partially confirmed
%                  1: confirmed

% path = truth.path;
% pos = truth.pos;
% seg_path = ones(2,length(all_pos));
% [Pos,idx1,idx2]=intersect(all_pos,pos);
% st=1;
% for i = 2:length(idx1)
%   ii = find(model.a_trans(1,:)==(path(st));
%   ii2 = find(model.a_trans(2,ii)==(path(i));
%   seg_path(st:idx(i)) = model.a_trans(3,ii(ii2));
%   ii = find(all_pos>= Pos(i)- model.loss.zero_range(path(st),path(i))...
%              &all_pos<= Pos(i)+ model.loss.zero_range(path(st),path(i));
%   seg_path(2,idx1(i)) = 1; 
% end

seg = truth.segments;
%seg_path = ones(2,length(all_pos));
seg_path = zeros(2,max(all_pos));

for i=1:size(seg,1)
  seg_path(1, seg(i,1):seg(i,2)-1) = seg(i,3) ;
end ;

% confirmed genes are worth more than unconfirmed ones
seg_path(2,:) = model.loss.unconfirmed ;

gene_idx = setdiff(unique(truth.segments(:,4)),0);
for j=1:length(gene_idx)
  idx = find(truth.segments(:,4)==gene_idx(j));
  genes{j} = truth.segments(idx,1:3) ;
end

%assert(length(genes)==length(gene_confirmed)) ;
for i=1:length(genes)
 % if gene_confirmed(i)==1%confirmed
    start = genes{i}(1,1) ;
    stop  = genes{i}(end,2) ;
    seg_path(2, max(1,start-model.loss.confirmed_tol) : min(size(seg_path,2), stop+model.loss.confirmed_tol) ) = model.loss.confirmed ;
  %elseif gene_confirmed(i)==0%partially confirmed 
  %  start = genes{i}(1,1) ;
  %  stop  = genes{i}(end,2) ;
  %  seg_path(2, max(1,start-model.loss.confirmed_tol) : min(size(seg_path,2), stop+model.loss.confirmed_tol) ) = (model.loss.confirmed+model.loss.unconfirmed)/2 ;
  %end ;
end ;

pos = truth.pos;
[Pos,idx1,idx2]=intersect(all_pos, pos);
assert(length(Pos)==length(pos)) ;
signal_names = fieldnames(model.signals) ;
for i = 1:length(idx1)
  signal = model.states(truth.path(idx2(i))).signal ;
  if length(signal)~=1, continue ; end ;
  found = 0 ;
  for j=1:model.cnt_signals,
    if getfield(model.signals, signal_names{j})==signal
      signal = j ;
      found = 1 ;
      break ;
    end ;
  end ;
  assert(found~=0) ;
  tol = getfield(model.loss.zero_range, signal_names{signal}) ;
  if tol==0, continue ; end ;
  idx = max(1,(Pos(i) - tol)):min((Pos(i)+ tol),size(seg_path,2)) ;
  seg_path(2,idx) = 0 ; 
end

seg_path = seg_path(:, all_pos) ;

truth.seg_path = seg_path ;
