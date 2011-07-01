function QP = solve_qp_iter(QP)

for j = [1 5 20 50 100]
	fprintf('solving QP using %i%% percent of the inactive constraints\n', j);
	[QP frac_solved] = solve_qp(QP, [], j) ;
	if frac_solved>0.99
		return
	end
end
