function my_hist(arr, col, num_bins)

if nargin<2
	col = 'b';
end
if nargin<3
	num_bins = 100;	
end

min_ = min(arr);
max_ = max(arr);

step = (max_-min_)/num_bins;
bins = min_:step:max_;

hi = histc(arr, bins);

hi = hi/sum(hi);

plot(bins, hi, col);
set(gca, 'XTickLabel', exp(bins(1:10:100)))
