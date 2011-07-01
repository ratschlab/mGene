function [XT,LT,POS,SC] = subsample(XT, LT, POS, SC, subsample_neg, subsample_pos)
% [XT,LT,POS,SC] = subsample(XT, LT, POS, SC, subsample_neg, subsample_pos)

if subsample_neg<1 & subsample_neg>0 & subsample_pos<=1 & subsample_pos>0 % & ~isempty(LT) 

  fprintf('subsampling %i %% of negative examples \n',subsample_neg*100);
  fprintf('subsampling %i %% of positive examples \n',subsample_pos*100);
  posidx = find(LT==1);
  negidx = find(LT==-1);
  rand('seed', pi) ;
  idx2 = randperm(length(negidx));
  idx3 = randperm(length(posidx));
  num_neg = round(subsample_neg*sum(LT==-1));
  num_pos = round(subsample_pos*sum(LT==1));
  
  XT = [XT(:,posidx(idx3(1:num_pos))) XT(:,negidx(idx2(1:num_neg)))] ;
  
  if ~isempty(SC)
    SC = [SC(:,posidx(idx3(1:num_pos))) SC(:,negidx(idx2(1:num_neg)))] ;
  end
  LT = [LT(posidx(idx3(1:num_pos))) LT(negidx(idx2(1:num_neg)))] ; 
  POS.pos = [POS.pos(posidx(idx3(1:num_pos))) POS.pos(negidx(idx2(1:num_neg)))] ; 
  POS.region_id = [POS.region_id(posidx(idx3(1:num_pos))) POS.region_id(negidx(idx2(1:num_neg)))] ; 
  fprintf('%i examples remaining (%.2f %% positives)\n',length(LT),sum(LT==1)/length(LT)*100);

end
