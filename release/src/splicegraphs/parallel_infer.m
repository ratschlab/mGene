addpath(sprintf('%s/splicegraphs', origpath)) ;
addpath /fml/ag-raetsch/share/software/matlab_tools/utils/
addpath ~/svn/tools/rproc/


num_jobs = 37;
chunksize = round(length(genes)/num_jobs);
MEMREQ = 3000;     % approximate memory requirement in MB
TIMEREQ = 500;     % approximate time requirement in minutes


joblist=rproc_empty(0) ;
options={} ;
options.identifier = 'inf';
options.start_dir = [getenv('HOME') '/svn/projects/splicing/splicegraphs'] ;



gene_subset = [1:chunksize];
fprintf(1,'Submitting infer(1:%d)\n',chunksize);
joblist = [joblist, rproc('infer',genes(gene_subset),MEMREQ,options,TIMEREQ) ];

%submit a job for each chunk
for ix=2:num_jobs-1
  gene_subset = [(ix-1)*chunksize+1:ix*chunksize];
  fprintf(1,'Submitting infer(%d:%d)\n',(ix-1)*chunksize+1,ix*chunksize);
  joblist = [joblist, rproc('infer',genes(gene_subset),MEMREQ,options,TIMEREQ) ];
end;


gene_subset = [(num_jobs-1)*chunksize+1:length(genes)];
fprintf(1,'Submitting infer(%d:%d)\n',(num_jobs-1)*chunksize+1,length(genes));
joblist = [joblist, rproc('infer',genes(gene_subset),MEMREQ,options,TIMEREQ) ];

% wait till all are finished
[joblist, num_crashed] = rproc_wait(joblist);
if num_crashed > 0
  fprintf(1,'Error: jobs crashed\n');
end


genes = [] ;
for ix=1:num_jobs
  fprintf(1,'Collecting %d\n',ix);
  genes = [genes,rproc_result(joblist(ix))];
end;

%fprintf(1,'Saving results\n');
%save /fml/ag-raetsch/home/ong/MatlabCode/SpliceGraph/Data/elegans.genes.p_infer.mat genes;

rproc_cleanup(joblist);
clear joblist;

