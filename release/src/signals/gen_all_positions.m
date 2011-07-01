function  jobinfo = gen_all_positions(fn_genome_config,fn_pos,Signal,fn_training_blocks,fn_test_blocks, num_splits,run_locally,RPROC) 
  
jobinfo = rproc_empty(0);
if ~run_locally,  
  jobinfo = rproc_empty(0) ;
  [sg_path, envstr] = shogun_settings();
  options = RPROC.options;
  options.start_dir = get_base_dir();
  options.addpaths = {sg_path, fileparts(which('gen_all_positions_rproc'))};
  options.envstr = envstr;
  options.ncpus = 1;
  options.hard_time_limit = 10;%kill job after 10 minutes
  % options.identifier = sprintf('%spos%s_%s',upper(PAR.organism.name(1:4)),PAR.Signal_name,PAR.method.name) ;
  % options.identifier = sprintf('pos%s_%s',PAR.Signal_name,PAR.method.name) ;
end

strands = '+-';
genome_info = init_genome(fn_genome_config);

% keyboard

blocks = [];
if num_splits~=1 
  % load(FN.output.fn_regions, 'regions')
  % load(FN.output.fn_split,'split')
  if ~isempty(fn_training_blocks) && fexist(fn_training_blocks)

    % load training splits
    TR = [] ;
    TR.blocks = load_struct(fn_training_blocks, 'blocks');
    tmp = load(fn_training_blocks, 'split') ;
    TR.split = tmp.split ;
    if ~isempty(TR.split)& ~isempty(TR.blocks)
      has_split = isfield(TR.blocks,'split')&&~isempty(TR.blocks(1).split);
      splits = fieldnames(TR.split);
      for f=1:length(splits)
        regids = getfield(TR.split, splits{f});
        [temp i1 i2] = intersect(regids, [TR.blocks.id]);
        if ~has_split
          for ii2=i2, 
            [TR.blocks(ii2).split] = deal(f);
          end 
          %[TR.blocks(i2).split] = deal(f);
        else
          assert(all([TR.blocks(i2).split] == f))      
        end
      end
      clear temp i1 i2 splits regids 
      blocks = [blocks TR.blocks];
      clear TR
    end
  end
  if ~isempty(fn_test_blocks) && fexist(fn_test_blocks)
    TE=[] ;
    TE.blocks = load_struct(fn_test_blocks, 'blocks') ;
    tmp = load(fn_test_blocks, 'split');
    TE.split = tmp.split ;
    has_split = isfield(TE.blocks,'split');  
    splits = fields(TE.split);
    for f=1:length(splits)
      regids = getfield(TE.split, splits{f});
      [temp i1 i2] = intersect(regids, [TE.blocks.id]);
      if ~has_split
        for ii2=i2, 
          [TE.blocks(ii2).split] = deal(f);
        end 
        %[TE.blocks(i2).split] = deal(f);
      else
        assert(all([TE.blocks(i2).split] == f))      
      end
    end
    if isfield(TE.blocks,'right_dist')
      TE.blocks = rmfield(TE.blocks,'right_dist');
    end
    if isfield(TE.blocks,'right_dist')
      TE.blocks = rmfield(TE.blocks,'right_dist');
    end
    blocks = [blocks TE.blocks];
    clear TE;
  end
end

%[temp, num_files] = unix(sprintf('ls %scontig_*.pos | wc -l',fn_pos));
%d=dir(sprintf('%scontig_*.pos',fn_pos));
%num_files = length(d) ;
%[temp, num_files1] = unix(sprintf('ls %scontig_*.svm | wc -l',fn_pos));
%d=dir(sprintf('%scontig_*.svm',fn_pos));
%num_files1=length(d) ;

%----------
% generate positions
%-----------
trivial_regions = init_regions(fn_genome_config);

% comp_region  = ones(1,length(trivial_regions));
num_jobs=0;


for t=1:length(trivial_regions)
  if fexist(sprintf('%scontig_%i%s.pos',fn_pos,trivial_regions(t).chr_num,trivial_regions(t).strand))&...
        fexist(sprintf('%scontig_%i%s.svm',fn_pos,trivial_regions(t).chr_num,trivial_regions(t).strand))
    % comp_region(t)  = 0;
    continue
  end
  num_jobs = num_jobs+1;
  
  P.region = trivial_regions(t); 
  P.blocks = blocks;
  P.Signal = Signal;
  P.fn_pos = fn_pos;
  P.num_splits = num_splits;

  if ~run_locally
    [mem_req, time_req, opts] = rproc_memtime_policy('gen_all_positions_rproc',  0, options) ;
    jobinfo(num_jobs) = rproc_create('gen_all_positions_rproc', P, mem_req, opts, time_req) ;
  else
    gen_all_positions_rproc(P);
  end
end

return
