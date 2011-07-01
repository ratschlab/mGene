function vec=subsample_policy(size_vec, sensor_name)
% input: 	size_vec: 	generic vector of sample sizes with dimension depending on sensor_name
% 		sensor_name: 	possible values are 'acc', 'don', 'tis', 'cdsStop', 'tss', 'cleave', 'polya', 
%				'cds_exon', 'frame0', 'intergenic', 'intron', 'polya_tail',  'utr3exon', 'utr5exon' 


num_splits = 5;

if strcmp(sensor_name, 'acc')||strcmp(sensor_name, 'don')
  	%maxval = floor(5e5/num_splits);
  	maxval = floor(5e6/num_splits);
	%vec = size_vec; 
	%return
elseif strcmp(sensor_name, 'tis')||strcmp(sensor_name, 'cdsStop')||...
	strcmp(sensor_name, 'tss')||strcmp(sensor_name, 'cleave')||...
	strcmp(sensor_name, 'polya')
  	maxval = floor(1e5/num_splits);
	%vec = size_vec;
	%return
elseif strcmp(sensor_name, 'cds_exon')||strcmp(sensor_name, 'frame0')||...
	strcmp(sensor_name, 'intergenic')||strcmp(sensor_name, 'intron')||...
	strcmp(sensor_name, 'polya_tail')||strcmp(sensor_name, 'utr3exon')||...
	strcmp(sensor_name, 'utr5exon')
  maxval = floor(1e5/num_splits);
end

if strcmp(sensor_name, 'acc')||strcmp(sensor_name, 'don')||...
   	strcmp(sensor_name, 'tis')||strcmp(sensor_name, 'cdsStop')||...
   	strcmp(sensor_name, 'tss')||strcmp(sensor_name, 'cleave')||...
   	strcmp(sensor_name, 'polya')||...
	strcmp(sensor_name, 'cds_exon')||strcmp(sensor_name, 'frame0')||...
	strcmp(sensor_name, 'intergenic')||strcmp(sensor_name, 'intron')||...
	strcmp(sensor_name, 'polya_tail')||strcmp(sensor_name, 'utr3exon')||...
	strcmp(sensor_name, 'utr5exon')

  assert(length(size_vec)==2)
  num_pos = size_vec(1);
  num_neg = size_vec(2);

  if num_pos+num_neg<maxval
    vec = [num_pos num_neg];
  elseif num_pos ==0
    vec = [num_pos min(maxval, num_neg)];
  elseif num_pos*11 < maxval
    vec = [num_pos min(10*num_pos, num_neg)];
  elseif maxval/num_pos>3
    vec = [num_pos min(maxval-num_pos, num_neg)];
  else
    vec = [floor(maxval/3) min(num_neg, floor((2*maxval)/3))];
  end
elseif strcmp(sensor_name, 'gene')
  assert(length(size_vec)==1)
  num_blocks = size_vec;
  maxval = 5000;
  minval = 100;

  if num_blocks>maxval
    num_train = floor(maxval/2);
    num_test = floor(maxval/2);
    %num_test = min(max_val, (num_blocks-num_train);
  elseif 0%num_blocks>minval
    num_train = floor((3*num_blocks)/4);
    num_test = floor(num_blocks/4);
  else
    num_train = num_blocks;
    num_test = 0;
  end
  vec = [num_train num_test];
  assert(num_train+num_test<=num_blocks)
else
  error('could not handle input argument: %s ', sensor_name);
end 
