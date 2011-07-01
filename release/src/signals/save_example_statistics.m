function stat=save_example_statistics(Signal,fn_filter_settings,fn_examples,fn_example_statistics,num_splits)
% stat = save_example_statistics(Signal,fn_filter_settings,fn_examples,fn_example_statistics,num_splits)


pos = 0;
neg = 0;
for j=1:num_splits
  f = gen_example_filename(Signal,fn_filter_settings,fn_examples,1,j) ;
  l=load(f, 'LT');
  pos = pos + sum(l.LT==1);
  neg = neg + sum(l.LT==-1);
  stat.(sprintf('split_%i_num_pos',j))=sum(l.LT==1);
  stat.(sprintf('split_%i_num_neg',j))=sum(l.LT==-1);
  stat.(sprintf('split_%i_frac_pos',j))=sum(l.LT==1)/sum(l.LT==-1)*100;
end

stat.all_num_pos = pos;
stat.all_num_neg = neg;
stat.all_frac_pos=pos/neg*100;

assert(pos>0);
assert(neg>0);


fn = fn_example_statistics(1:end-4);
num = 1;
while fexist(sprintf('%s_%i.mat',fn,num))
  num=num+1;
end

save(sprintf('%s_%i.mat',fn,num), 'stat');
