function include_preds_all(origpath,organismlistfile)
% A wrapper script to include the predicted alternative splicing events
% into the splice graph

if nargin < 2
  organismlistfile = 'orglist.txt';
end
if nargin < 1
  origpath = '~/svn/projects/splicing';
end

addpath ~/svn/tools/rproc
addpath ~/svn/tools/utils

orglist = textread(organismlistfile,'%s')
numorgs = length(orglist);


% setup rproc
options={} ;
options.identifier = 'PA';
options.start_dir = [getenv('HOME') '/svn/projects/splicing/splicegraphs'] ;
MEMREQ = 5000;

joblist=rproc_empty(0) ;
% do the includes for each fdr
for ixr = [5 10 25 50]
  for ixo = 1:numorgs
    data.ginfo = sprintf('/fml/ag-raetsch/share/projects/altsplicedata/%s/genome.config',orglist{ixo});
    data.origpath = origpath;
    data.fd_rate = ixr;
    if strcmp(orglist{ixo},'human')
      joblist(end+1) = rproc('include_preds',data,15000,options,30);
    else
      joblist(end+1) = rproc('include_preds',data,MEMREQ,options,30);
    end
    pause(1);
  end
end

[joblist, num_crashed] = rproc_wait(joblist,60,1,0);
if num_crashed > 0
  fprintf(1,'Error: jobs crashed\n');
  [joblist, num_crashed] = rproc_wait(joblist,60,1,0);
  %keyboard;
end

if num_crashed == 0
  rproc_cleanup(joblist);
  clear joblist;
end
