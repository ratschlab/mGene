function x=save_label_files_starter(D)
% save_label_files_starter(D)

% use the rproc.addpath way to add paths in client processes!!
%paths;

region = D.t_region;
region = load_sequence(region) ;

P = D.P;

if ~P.label_source.from_fn_candsites,

  genes = load_struct(P.fn_genes_signals,'genes') ;
  
  save_label_files(region, P, genes) ;
else
  save_label_files(region, P) ;
end ;

% dummy return argument (rproc requirement)
x=0;
