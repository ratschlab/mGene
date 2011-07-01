function [total_num_pos, num_pos] = retrieve_num_pos(regions,base_dir)
% num_pos = retrieve_num_pos(regions,base_dir)
%

num_pos = zeros(1,length(regions)) ;
for r=1:length(regions)
  f_name = sprintf('%scontig_%i%s',base_dir,regions(r).chr_num,regions(r).strand);
  d=dir([f_name '.pos']) ;
  if d(1).bytes==0,
    warning('empty pos file') ;
    num_pos(r) = 0 ;
  else
    num_pos(r) = get_num_pos(f_name, [regions(r).start; regions(r).stop]) ;
  end ;
end
total_num_pos = sum(num_pos) ;
