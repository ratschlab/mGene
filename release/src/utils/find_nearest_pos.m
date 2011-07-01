function pos=find_nearest_pos(old_pos, all_pos, lb, ub)

	cands = all_pos(find(all_pos<=ub&all_pos>=lb));
	if isempty(cands), 
		pos=[] ; 
		return ; 
	end ;
	[mm idx] = min(abs(cands-old_pos));
	pos = cands(idx);
return
