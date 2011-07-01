function sg_set_plif(penalty_array)
% sg_set_plif(penalty_array)

for j=1:length(penalty_array)
  all_ids(j)		= penalty_array{j}.id;
  all_names{j} 		= penalty_array{j}.name;
  all_limits(:,j) 	= penalty_array{j}.limits';
  all_penalties(:,j) 	= penalty_array{j}.penalties';
  if isempty(penalty_array{j}.transform)
    all_transform{j} 	= 'linear';
  else
    all_transform{j} 	= penalty_array{j}.transform;
  end
  all_min_values(j) 	= penalty_array{j}.min_value;
  all_max_values(j) 	= penalty_array{j}.max_value;
  all_use_cache(j) 	= penalty_array{j}.use_cache;
  all_use_svm(j) 	= penalty_array{j}.use_svm;
  all_do_calc(j) 	= 1;%penalty_array{j}.do_calc;
end

sg('set_plif_struct',int32(all_ids)-1,all_names, all_limits, all_penalties, all_transform, all_min_values, all_max_values, int32(all_use_cache), int32(all_use_svm), int32(all_do_calc)); 
[ret_ids,ret_names, ret_limits, ret_penalties, ret_transform, ret_min_values, ret_max_values, ret_use_cache, ret_use_svm, ret_do_calc]= sg('get_plif_struct');

[engine, env]=determine_engine() ;

assert(isequal(ret_do_calc, all_do_calc))
assert(isequal(ret_use_svm, all_use_svm))
assert(isequal(ret_use_cache, all_use_cache))
assert(isequal(ret_max_values, all_max_values))
assert(isequal(ret_min_values, all_min_values))
if ~isequal(engine, 'octave'),
  assert(isequal(ret_transform', all_transform))
end ;
assert(isequal(ret_penalties, all_penalties))
assert(isequal(ret_limits, all_limits))
if ~isequal(engine, 'octave'),
  assert(isequal(ret_names', all_names))
end ;
assert(isequal(ret_ids+1, all_ids))

