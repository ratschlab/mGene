function PAR=masterscript_data(PAR,first)
  
%PAR=masterscript_data(PAR,first)

if nargin<2
  first =1 ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREPARE COMPLETE GENE STRUCTURE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PAR.tasks.prepare_labels
  if first==1 
    fprintf('start preparing gene structure for organism: %s\n\n',PAR.organism.name)
    PAR.add_signals.tss = 1; 
    PAR.add_signals.cleave = 1 ; 
    PAR.add_signals.polya = 1;
    PAR.add_signals.tis = 1;
    PAR.add_signals.cdsStop = 1;
    PAR.add_signals.transacc = 0; 
  elseif first==2
    fprintf('adding transacc sites to gene structure for organism: %s\n\n',PAR.organism.name)
    PAR.add_signals.tss = 0; 
    PAR.add_signals.cleave = 0 ; 
    PAR.add_signals.polya = 0;
    PAR.add_signals.tis = 0;
    PAR.add_signals.cdsStop = 0;
    PAR.add_signals.transacc = 1; 
  end
  
  [genes_sensors,genes_merged] = prepare_complete_genes(PAR);
  if first==2
    PAR.add_signals.tss = 1; 
    PAR.add_signals.cleave = 1 ; 
    PAR.add_signals.polya = 1;
    PAR.add_signals.tis = 1;
    PAR.add_signals.cdsStop = 1;
    PAR.add_signals.transacc = 1; 
  end
  genes = genes_sensors;
  if fexist(PAR.FN.output.fn_genes_signals) || fexist([PAR.FN.output.fn_genes_signals '.mat']),
    warning(sprintf('do not save file as it already exists: %s', PAR.FN.output.fn_genes_signals)) ;
    L=load_struct(PAR.FN.output.fn_genes_signals, 'genes') ;
    assert(isequal(L, genes)) ;
    clear L
  else
    save_struct(PAR.FN.output.fn_genes_signals,genes,'genes')
    %save(PAR.FN.output.fn_genes_signals,'-append','PAR')
    save_append(PAR.FN.output.fn_genes_signals, 1, 'PAR', PAR)
  end ;
  genes = genes_merged;
  if fexist(PAR.FN.output.fn_genes_merged) || fexist([PAR.FN.output.fn_genes_merged '.mat']),
    warning(sprintf('do not save file as it already exists: %s', PAR.FN.output.fn_genes_signals)) ;
    L=load_struct(PAR.FN.output.fn_genes_merged, 'genes') ;
    assert(isequal(L, genes)) ;
    clear L
  else
    save_struct(PAR.FN.output.fn_genes_merged,genes,'genes')
    %save(PAR.FN.output.fn_genes_signals,'-append','PAR')
    save_append(PAR.FN.output.fn_genes_signals, 1, 'PAR', PAR)
  end 
  clear genes_merged, genes_sensors
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREPARE REGIONS AND BLOCKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% keyboard
if first==1 &...
         (~fexist(PAR.FN.output.fn_trivial_regions) |...
          ~fexist(PAR.FN.output.fn_blocks)|...
          ~fexist(PAR.FN.output.fn_test_blocks))
  [trivial_regions,blocks,training,test] = prepare_regions(PAR);
  if ~fexist(PAR.FN.output.fn_trivial_regions)| ~fexist(PAR.FN.output.fn_blocks)  
    save(PAR.FN.output.fn_trivial_regions,'trivial_regions','PAR');
    save(PAR.FN.output.fn_blocks, 'blocks','PAR');
  end
  if 1 % ~isempty(training.blocks)&~fexist(PAR.FN.output.fn_training_blocks)
    %save(PAR.FN.output.fn_training_blocks,'-struct','training');
    save_append(PAR.FN.output.fn_training_blocks, 0, training);
    %save(PAR.FN.output.fn_training_blocks,'-append' ,'PAR');
    save_append(PAR.FN.output.fn_training_blocks, 1, 'PAR', PAR);
    PAR.regions.offset=PAR.regions.offset_train;
    save_splits2bed(training.blocks,training.split,PAR.FN.output.fn_training_blocks(1:end-4),PAR);
  end
  if 1 % ~isempty(test.blocks)&~fexist(PAR.FN.output.fn_test_blocks)
    %save(PAR.FN.output.fn_test_blocks,'-struct','test');
    save_append(PAR.FN.output.fn_test_blocks, 0, test);
    %save(PAR.FN.output.fn_test_blocks,'-append','PAR');
    save_append(PAR.FN.output.fn_test_blocks, 1, 'PAR', PAR);
    PAR.regions.offset=PAR.regions.offset_test;
    save_splits2bed(test.blocks,test.split,PAR.FN.output.fn_test_blocks(1:end-4),PAR);
  end
end
