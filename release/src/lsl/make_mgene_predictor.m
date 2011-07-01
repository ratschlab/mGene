function make_mgene_predictor(predictor_fn, predictor_dir,...
			tss_fn, tss_dir, tis_fn, tis_dir, acc_fn, acc_dir, don_fn, don_dir, cdsStop_fn, cdsStop_dir, cleave_fn, cleave_dir,...
			intergenic_fn, intergenic_dir, utr5exon_fn, utr5exon_dir, cds_exon_fn, cds_exon_dir, intron_fn, intron_dir, utr3exon_fn, utr3exon_dir,...
			fn_out, output_dir)

fd=fopen(fn_out, 'w+') ;
[ret, timedate] = unix('date') ;
fprintf(fd, 'mGene Predictor created on %s\n\n', timedate)
fclose(fd) ;

%signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
train_dirs.tss = tss_dir;
train_files.tss = tss_fn;

train_dirs.tis = tis_dir;
train_files.tis = tis_fn;

train_dirs.acc = acc_dir;
train_files.acc = acc_fn;

train_dirs.don = don_dir;
train_files.don = don_fn;

train_dirs.cdsStop = cdsStop_dir;
train_files.cdsStop = cdsStop_fn;

train_dirs.cleave = cleave_dir;
train_files.cleave = cleave_fn;


%contents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
train_dirs.intergenic = intergenic_dir;
train_files.intergenic = intergenic_fn;

train_dirs.utr5exon = utr5exon_dir;
train_files.utr5exon = utr5exon_fn;

train_dirs.cds_exon = cds_exon_dir;
train_files.cds_exon = cds_exon_fn;

train_dirs.intron = intron_dir;
train_files.intron = intron_fn;

train_dirs.utr3exon = utr3exon_dir;
train_files.utr3exon = utr3exon_fn;


%genes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
train_dirs.gene_predictor = predictor_dir;
train_files.gene_predictor = predictor_fn;

fn={'gene_predictor', 'tss', 'tis', 'acc', 'don', 'cdsStop', 'cleave', ...
   'intergenic', 'utr5exon', 'cds_exon', 'intron', 'utr3exon'} ;
for i=1:length(fn)
  system(sprintf('cat %s >> %s && echo >> %s && echo ------ >> %s && echo >> %s ', train_files.(fn{i}), ...
                 fn_out, fn_out, fn_out, fn_out)) ;
end ;

nice_mkdir(output_dir) ;

save(sprintf('%s/mgene_predictor.mat', output_dir), 'train_dirs', ...
     'train_files', '-v7');


%function fuse_mGene_predictor(	fn_fused_mGene_predictor, fn_gene_predictor,...
%				tss_dir, tis_dir, acc_dir, don_dir, stop_dir, cleave_dir,...
%				exon_dir, intron_dir, utr_dir, intergenic_dir)


%% collect signal classifiers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%signal_dirs = {tss_dir, tis_dir, acc_dir, don_dir, stop_dir, cleave_dir};
%
%for k=1:length(signal_dirs)
%  l = load([signal_dirs{k} '/PAR.mat'], 'PAR');
%  signal = l.PAR.Signal_name;
%  for j=1:5
%    fn_SVM = sprintf('%sbest_partition=%i.mat', l.PAR.FN.output_sig.(signal).fn_SVMs, j);
%    signal_classifiers.(signal) = load(fn_SVM);
%  end
%end
%
%clear l fn_SVM
%% collect content classifiers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%content_dirs = {exon_dir, intron_dir, utr_dir, intergenic_dir}
%
%for k=1:length(content_dirs)
%  if fexist([content_dirs{k} '/PAR.mat'])
%    l = load([content_dirs{k} '/PAR.mat'], 'PAR');
%    content = l.PAR.Content_name;
%    for j=1:5
%      fn_SVM = sprintf('%sbest_partition=%i.mat', l.PAR.FN.output_cont.(content).fn_SVMs, j);
%      if fexist(fn_SVM)
%        content_classifiers.(content) = load(fn_SVM);
%      end
%    end
%  end
%end
%
%% get gene structure predictor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%gene_predictor = struct;
%if fexist(fn_gene_predictor)
% gene_predictor = load(fn_gene_predictor);
%else
%  warning('no gene predictor found')
%end
%
%save(fn_fused_mGene_predictor, 'content_classifiers', 'signal_classifiers', 'gene_predictor')
%
%
%
%
