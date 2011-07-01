function sp = sparse_features(features)

assert(all(all(all(features(~isinf(features))>=0))))
features = features+1;
features(isinf(features))=0;

sp{1} = sparse(features(:,:,1));
sp{2} = sparse(features(:,:,2));

