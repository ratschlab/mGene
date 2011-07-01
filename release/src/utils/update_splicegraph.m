function genes = update_splicegraph(genes, run_locally)

%addpath ~/svn/projects/genefinding
%paths
if nargin==1
  run_locally=1;
end

if run_locally
  genes = rmfield(genes, 'splicegraph');
  genes = build_splice_graph_caller(genes);
  genes = infer_splice_graph_caller(genes);
  genes = alt_const(genes);
  %genes = correct_splicegraph(genes);
else
  for j=1:length(genes)
    dir_ = which('update_splicegraph');
    dir_(find(dir_=='/',1,'last')+1:end)=[];
    options.start_dir = dir_;
    jobinfo(j)=rproc_create('update_splicegraph',genes(j),500,options,60);
  end
  jobinfo = rproc_submit_batch(jobinfo,1000);
  rproc_wait(jobinfo, 60,1,0);
  for j = 1:length(genes)
    l = load(jobinfo(j).result_fname)
    genes(j)=l.retval1;
    clear l
  end
end

