function FN=init_FN(signal_names,content_names)
 
  
if nargin<1
  signal_names = {'acc','don','tis','cdsStop','tss','cleave','polya','transacc'};
end
if nargin<2
  content_names = {'intergenic','utr5exon','cds_exon','utr3exon','polya_tail','intron','frame0'};
end
 
%--------------------------------------------------------------------------
% Directories
% -------------------------------------------------------------------------------------

FN.genome_dir= '/fml/ag-raetsch/share/databases/genomes/';
FN.input_dir= '/fml/ag-raetsch/share/projects/genefinding/';
FN.output_dir = '/fml/ag-raetsch/share/projects/genefinding/';

FN.dir_name = '';
FN.exp_name = '';
FN.lsl_exp_name = '';
%--------------------------------------------------------------------------
% INPUT FILE NAMES
% -------------------------------------------------------------------------------------

FN.genome.fn_confirmed_sequences = 'confirmed_sequences.mat';
FN.genome.fn_genome_config_orig  = 'genome_uncondensed.config';
FN.genome.fn_genome_config  = 'genome.config';
FN.genome.fn_map = 'condense_map.txt';
FN.genome.fn_genes_anno = 'annotated_genes.mat';
FN.genome.fn_genes_merged = '';
FN.genome.fn_genes_conf = 'confirmed_sequences.mat';
FN.genome.fn_name_to_contig_mapping = 'sequences/contig_orig_new_mapping';


FN.input.fn_genome_config = 'genome.config';
FN.input.fn_trainings_regions = 'training_regions.mat';
FN.input.fn_test_regions = 'test_regions.mat';
FN.input.fn_genes_anno_orig ='annotated_genes.mat';
FN.input.fn_genes_conf_orig = 'confirmed_sequences.mat';
FN.input.fn_ESTs_matches = {'est_matches.mat'};

%--------------------------------------------------------------------------
% OUTPUT FILE NAMES
% -------------------------------------------------------------------------------------

FN.output.fn_genes_anno ='anno_genes.mat';
FN.output.fn_genes_conf = 'conf_genes.mat';
FN.output.fn_genes_merged = 'genes_all.mat';
FN.output.fn_genes_signals = 'genes_signals.mat';
FN.output.fn_genes_contents = 'annotated_genes.mat';
FN.output.fn_genes_lsl = 'annotated_genes.mat';
FN.output.fn_ESTs_matches = [];
FN.output.fn_trivial_regions = 'trivial_regions.mat';
FN.output.fn_blocks = 'blocks.mat';
FN.output.fn_training_blocks ='training_blocks.mat';
FN.output.fn_test_blocks = 'test_blocks.mat';

%--------------------------------------------------------------------------
% SIGNAL FILE NAMES
% -------------------------------------------------------------------------------------

for j=1:length(signal_names)
  sig_name = signal_names{j};
  FN.input_sig.(sig_name).fn_candsites = 'cands/' ;
  FN.input_sig.(sig_name).fn_examples = 'examples';
  FN.input_sig.(sig_name).fn_pos = 'pos/';
  FN.input_sig.(sig_name).fn_filter_settings= 'filter_settings.mat';
  FN.input_sig.(sig_name).fn_example_statistics = 'example_statistics.mat';

  FN.output_sig.(sig_name).fn_models='MS/models';
  FN.output_sig.(sig_name).fn_SVMs='MS/';
  FN.output_sig.(sig_name).fn_pred='pred/';
  FN.output_sig.(sig_name).fn_motifs='MS/motifs';
end


if any(strcmp(signal_names,'tss'))
  FN.input_sig.tss.fn_est_cluster = 'tss_est_clusters.mat';
  FN.input_sig.tss.fn_fulllength  = 'tss_fulllength.mat';
  FN.input_sig.tss.fn_db = 'tss_db.mat' ;
end

if any(strcmp(signal_names,'cleave'))
  FN.input_sig.cleave.fn_est_cluster = 'cleave_est_clusters.mat';
  FN.input_sig.cleave.fn_fulllength = 'cleave_fulllength.mat';
  FN.input_sig.cleave.fn_db  = 'cleave_db.mat';
  FN.input_sig.cleave.fn_anno  = 'cleave_anno.mat';
end
if any(strcmp(signal_names,'polya'))
  FN.input_sig.polya.fn_polya = 'truesites_and_consensus.mat';
  FN.input_sig.polya.fn_log = 'log';
end
if any(strcmp(signal_names,'tis'))
  FN.input_sig.tis.fn_anno = 'tis_anno.mat';
end
if any(strcmp(signal_names,'cdsStop'))
  FN.input_sig.cdsStop.fn_anno = 'cdsStop_anno.mat';
end

if any(strcmp(signal_names,'transacc'))
  FN.input_sig.transacc.fn_anno_SL1 = 'transacc_anno_SL1' ;
  FN.input_sig.transacc.fn_anno_SL2 = 'transacc_anno_SL2' ;
  FN.input_sig.transacc.fn_pred = 'transacc_pred';
  FN.input_sig.transacc.fn_db = []; % only for C_elegans
end

%--------------------------------------------------------------------------
% CONTENT FILE NAMES
% -------------------------------------------------------------------------------------

for j=1:length(content_names)
  content_name = content_names{j};
  FN.input_cont.(content_name) = struct;
  FN.input_cont.(content_name).fn_candsites = 'cands/' ;
  FN.input_cont.(content_name).fn_examples = 'examples';
  FN.input_cont.(content_name).fn_pos = 'pos/';
  FN.input_cont.(content_name).fn_filter_settings= 'filter_settings.mat';
  FN.input_cont.(content_name).fn_example_statistics = 'example_statistics.mat';

  FN.output_cont.(content_name).fn_models = 'MS/models';
  FN.output_cont.(content_name).fn_SVMs = 'MS/';
  FN.output_cont.(content_name).fn_pred = 'pred/';
  FN.output_cont.(content_name).fn_motifs = 'MS/motifs';
end


%--------------------------------------------------------------------------
% LSL FILE NAMES
% -------------------------------------------------------------------------------------

FN.input_lsl.fn_test_blocks = 'data/test_blocks';
FN.input_lsl.fn_training_blocks =  'data/training_blocks.mat';
FN.input_lsl.fn_training_split =  'data/blocks_split.mat';

FN.input_lsl.fn_val_blocks =  'data/validation_blocks.mat';
FN.input_lsl.fn_blocks_all =  'data/blocks_all.mat';
FN.input_lsl.fn_boundary_model = 'data/boundary_model.mat';


FN.output_lsl.fn_init =  'train/init.mat';
FN.output_lsl.fn_train_block_preds = 'train/single_block_predictions/';
FN.output_lsl.fn_train = 'train/train';
FN.output_lsl.fn_pred = 'predictions/';
%FN.output_lsl.fn_blocks_pred = 'predictions/blocks_pred';
%FN.output_lsl.fn_pred_genome = 'genomewide_predictions/';
FN.output_lsl.fn_log = 'log_files/log_';


