function [l1, l2] = compute_overlap_lists(idx, num1, num2)
% idx is a Nx2 matrix indicating that region idx(i, 1) overlaps with region idx(i, 2)

l1 = cell(1, num1); 
l2 = cell(1, num2); 

for j = 1:size(idx, 1)
	l1{idx(j,1)} = [l1{idx(j,1)} idx(j,2)];
	l2{idx(j,2)} = [l2{idx(j,2)} idx(j,1)];
end
