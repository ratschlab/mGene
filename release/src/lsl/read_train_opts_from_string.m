function opts = read_train_opts_from_string(train_options, opts)

if ~isempty(train_options)
  train_option_items = separate(train_options, ';') ;
else
  train_option_items=[] ;
end ;
for i=1:length(train_option_items),
  if isempty(train_option_items{i}), 
    continue ;
  end ;
  items = separate(train_option_items{i},'=') ;
  if isequal(items{1}, 'block_design_merge_sep')
    opts.block_design_merge_sep = str2num(items{2}) ;
  elseif isequal(items{1}, 'maxNumIter')
    opts.maxNumIter = str2num(items{2}) ;
  elseif isequal(items{1}, 'maxNumBlock')
    if ~isequal(items{2},'-'),
      opts.maxNumBlock = str2num(items{2}) ;
      if opts.maxNumBlock==0,
        opts.maxNumBlock=inf ;
      end ;
    end ;
  elseif isequal(items{1}, 'C')
    opts.C = str2num(items{2}) ;
  elseif isequal(items{1}, 'long_trans_thresh')
    opts.long_trans_thresh = str2num(items{2});
  elseif isequal(items{1}, 'block_design') ;
    opts.block_design = items{2} ;
  elseif isequal(items{1}, 'block_design_merge_num') ;
    opts.block_design_merge_num = str2num(items{2}) ;
  elseif isequal(items{1}, 'use_loss_mask') ;
    opts.use_loss_mask = str2num(items{2}) ;
  elseif isequal(items{1}, 'cleave_offset') ;
    opts.cleave_offset = str2num(items{2}) ;
  elseif isequal(items{1}, 'scale_loss_mask') ;
    opts.scale_loss_mask = str2num(items{2}) ;
  elseif isequal(items{1}, 'block_design_regionlist')
    opts.block_design_regionlist = items{2};
  elseif isequal(items{1}, 'exon_map_file')
    opts.exon_map_file = items{2};
  elseif isequal(items{1}, 'subsample_reads')
    opts.subsample_reads = str2num(items{2});
  elseif isequal(items{1}, 'track_function')
    opts.track_functions{end+1} = items{2};
  elseif isequal(items{1}, 'track_file')
    opts.track_files{end+1} = items{2};
  elseif isequal(items{1}, 'track_param')
    opts.track_params{end+1} = items{2};
  elseif isequal(items{1}, 'track_name')
    opts.track_names{end+1} = items{2};
  elseif isequal(items{1}, 'track_monoton_function')
    opts.track_monoton_functions{end+1} = items{2};
  elseif isequal(items{1}, 'segment_feature_function')
    opts.segment_feature_functions{end+1} = items{2};
  elseif isequal(items{1}, 'segment_feature_file')
    opts.segment_feature_files{end+1} = items{2};
  elseif isequal(items{1}, 'segment_feature_param')
    opts.segment_feature_params{end+1} = items{2};
  elseif isequal(items{1}, 'segment_feature_name')
    opts.segment_feature_names{end+1} = items{2};
  elseif isequal(items{1}, 'segment_feature_monoton_function')
    opts.segment_feature_monoton_functions{end+1} = items{2};
  elseif isequal(items{1}, 'use_train_region')
    opts.use_train_region = items{2};
  elseif isequal(items{1}, 'train_region_file')
    opts.train_region_file = items{2};
  elseif isequal(items{1}, 'use_rna_seq_for_label_gen')
    opts.use_rna_seq_for_label_gen = str2num(items{2});
  elseif isempty(items{1})
  else
    error('unknown option: %s', train_option_items{i}) ;
  end ;
end ;

