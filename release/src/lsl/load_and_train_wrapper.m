function load_and_train_wrapper(fname)
% load_and_train_wrapper(fname)

load(fname, 'inventory') ;
for i=1:length(inventory),
	load(fname, inventory{i}) ;
end ;

if exist('Train_Data','var')
  x = load_and_train(P, Train_Data) ;
  inventory={'x'} ;
%  save(fname, 'inventory', 'P', 'Train_Data', 'x', '-V7') ;
  save(fname, 'inventory', 'x', '-V7') ;
else
  x = load_and_train(P) ;
  inventory={'x'} ;
%  save(fname, 'inventory', 'P', 'x', '-V7') ;
  save(fname, 'inventory', 'x', '-V7') ;
end ;

