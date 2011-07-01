
% Collect the statistics of alternative splicing prediction.
% How many of the alternative splicing events overlap?
% Two events overlap when they share a common exon.

info = '~/altsplicedata/test/genome.config';origpath='/fml/ag-raetsch/home/ong/splicing';

addpath(sprintf('%s/utils', origpath)) ;
addpath /fml/ag-raetsch/share/software/matlab_tools/utils


% --- initialise genome
datadir = 'prediction';

genome_info = init_genome(info) ;
fprintf('Loading genome at: %s\n', genome_info.basedir) ;
graphfilename=sprintf('%s/%s/cand_splicegraph.mat', genome_info.basedir, datadir) ;
load(graphfilename,'genes');
fprintf('loaded %s.\n', graphfilename) ;



% --- build summary statistics
for ix = 1:length(genes)
  genes(ix).maxevents = 0;
  if length(genes(ix).splicegraph) == 4
    isaltexon = zeros(4,size(genes(ix).splicegraph{1},2));
    for ixa = 1:length(genes(ix).alt_exon)
      isaltexon(1,genes(ix).alt_exon{ixa}) = isaltexon(1,genes(ix).alt_exon{ixa}) + 1;
    end
    for ixa = 1:length(genes(ix).alt_intron)
      isaltexon(2,genes(ix).alt_intron{ixa}) = isaltexon(2,genes(ix).alt_intron{ixa}) + 1;
    end
    for ixa = 1:length(genes(ix).alt_5prime)
      isaltexon(3,genes(ix).alt_5prime{ixa}) = isaltexon(3,genes(ix).alt_5prime{ixa}) + 1;
    end
    for ixa = 1:length(genes(ix).alt_3prime)
      isaltexon(4,genes(ix).alt_3prime{ixa}) = isaltexon(4,genes(ix).alt_3prime{ixa}) + 1;
    end
    genes(ix).isaltexon = isaltexon;
    genes(ix).maxevents = max(sum(isaltexon));
  end
end


% some code to view the splice graphs with predicted events
if false
  load ~/altsplicedata/test/prediction/cand_splicegraph.mat
  for ix = 1:length(genes)
    if length(genes(ix).splicegraph) == 4
      disp(ix);
      genes(ix).isaltexon
      viewsplicegraph_preds(genes(ix));
      keyboard;
    end
  end

  idx1 = find([genes(:).maxevents]==1);
  idx2 = find([genes(:).maxevents]==2);
  idx3 = find([genes(:).maxevents]==3);
  idx4 = find([genes(:).maxevents]>3);
  for ix = idx1
    disp(ix);
    genes(ix).isaltexon
    figure(1);
    viewsplicegraph_preds(genes(ix),1);
    figure(2);
    viewsplicegraph_preds(genes(ix));
    keyboard;
  end
  for ix = idx2
    disp(ix);
    genes(ix).isaltexon
    figure(1);
    viewsplicegraph_preds(genes(ix),1);
    figure(2);
    viewsplicegraph_preds(genes(ix));
    keyboard;
  end
  for ix = idx3
    disp(ix);
    genes(ix).isaltexon
    figure(1);
    viewsplicegraph_preds(genes(ix),1);
    figure(2);
    viewsplicegraph_preds(genes(ix));
    keyboard;
  end
  for ix = idx4
    disp(ix);
    genes(ix).isaltexon
    figure(1);
    viewsplicegraph_preds(genes(ix),1);
    figure(2);
    viewsplicegraph_preds(genes(ix));
    keyboard;
  end
end
