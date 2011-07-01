function count_alt_exons_all
% A wrapper script to count the fraction of alternative splicing
% Also includes predictions.

addpath ~/svn/tools/rproc
addpath ~/svn/tools/utils

origpath = '~/splicing';

% setup rproc
options={} ;
options.identifier = 'PA';
options.start_dir = [getenv('HOME') '/svn/projects/splicing/splicegraphs'] ;
MEMREQ = 5000;

joblist=rproc_empty(0) ;

include_preds_all(origpath,'orglist.txt');

data.organismlistfile = 'orglist.txt';
data.outfile = 'frac_alt_confirmed.txt';
data.summaryfile = 'frac_alt_confirmed_summary.txt';
data.sourcefile = 'confirmed';
data.fd_rate = 0;
joblist(end+1) = rproc('count_alt_exons',data,MEMREQ,options,30);
pause(1);

data.sourcefile = 'pred';
for fd_rate = [5, 10, 25, 50]
  data.outfile = sprintf('frac_alt_pred_%02d.txt',fd_rate);
  data.summaryfile = sprintf('frac_alt_pred_%02d_summary.txt',fd_rate);
  data.fd_rate = fd_rate;
  joblist(end+1) = rproc('count_alt_exons',data,MEMREQ,options,30);
  pause(1);
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

%include_preds_all(origpath,'mouse.txt');
%count_alt_exons('mouse.txt','frac_mouse_confirmed.txt','frac_mouse_confirmed_summary.txt','confirmed')
%for fd_rate = [5, 10, 25, 50]
%  count_alt_exons('mouse.txt',sprintf('frac_mouse_pred_%02d.txt',fd_rate),...
%		  sprintf('frac_mouse_pred_%02d_summary.txt',fd_rate),'pred',fd_rate);
%end


%include_preds_all(origpath,'human.txt');
%count_alt_exons('human.txt','frac_human_confirmed.txt','frac_human_confirmed_summary.txt','confirmed')
%for fd_rate = [5, 10, 25, 50]
%  count_alt_exons('human.txt',sprintf('frac_human_pred_%02d.txt',fd_rate),...
%		  sprintf('frac_human_pred_%02d_summary.txt',fd_rate),'pred',fd_rate);
%end

