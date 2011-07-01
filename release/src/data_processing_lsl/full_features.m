function features = full_features(sp)

assert(length(sp)==2)
assert(iscell(sp))
assert(isequal(size(sp{1}),size(sp{2})))

features = full(sp{1});
features(:,:,2) = full(sp{2});

features(features==0) = -inf;
features = features-1;
