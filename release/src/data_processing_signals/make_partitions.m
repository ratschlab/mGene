function SETs = make_partitions(SETs)
 

if SETs.num_partitions==1
  return
end
  
SETs.partitions = eye(SETs.num_partitions,SETs.num_partitions)+1;  
for i = 1:2:SETs.num_partitions-1
  SETs.partitions(i:i+1,i:i+1) = SETs.partitions(i:i+1,i:i+1)+1;
end
if mod(SETs.num_partitions,2)
  SETs.partitions(end,end-1:end) =  SETs.partitions(end,end-1:end)+1;
end
if SETs.num_splits~=SETs.num_partitions
  assert(mod(SETs.num_splits,SETs.num_partitions)==0)
  SETs.partitions = repmat(SETs.partitions,1,SETs.num_splits/SETs.num_partitions);
end

if ~isfield(SETs,'partition_train_on') 
  [tmp,idx] = unique(SETs.partitions==1,'rows');
  SETs.partition_train_on = zeros(1,SETs.num_partitions);
  SETs.partition_train_on(idx) = 1;
end



SETs.info.train = 1;
SETs.info.eval = 2;
SETs.info.plif = 2;
SETs.info.pred = 3;

for i=1:SETs.num_partitions
  assert(isequal(unique(SETs.partitions(i,:)),unique([SETs.info.train ...
                      SETs.info.eval SETs.info.plif SETs.info.pred ])))
end

assert(all(sum(SETs.partitions==SETs.info.pred)==1))
assert(length(unique(sum(SETs.partitions'==SETs.info.train)))==1)
assert(length(unique(sum(SETs.partitions'==SETs.info.eval)))==1)
assert(length(unique(sum(SETs.partitions'==SETs.info.plif)))==1)
assert(length(unique(sum(SETs.partitions'==SETs.info.pred)))==1)



