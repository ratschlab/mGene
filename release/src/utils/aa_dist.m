function distance = aa_dist(matches, tlen, qlen)
    distance = 1-(matches/max(tlen, qlen));
	assert(distance>=0)
	assert(distance<=1)
return

