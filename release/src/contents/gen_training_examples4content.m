function jobinfo = gen_training_examples4content(PAR)
% gen_labels_and_train(PAR)

%paths
  

Content_name = PAR.Content_name;
Content = PAR.Contents.(PAR.Content_name);
num_splits = PAR.SETs.num_splits;
label_source  = PAR.label_source.(Content_name);
info = PAR.info_genes.transcript;


FN = PAR.FN;
fn_genome_config = FN.input.fn_genome_config ;
fn_training_blocks = FN.output.fn_training_blocks;
fn_genes_contents = FN.output.fn_genes_contents ;

input_cont = FN.input_cont.(Content_name);
fn_candsites =  input_cont.fn_candsites;
fn_filter_settings = input_cont.fn_filter_settings;
fn_examples = input_cont.fn_examples;
fn_example_statistics = input_cont.fn_example_statistics;

run_locally = PAR.tasks.contents.run_locally.pos ;

if ~run_locally
  RPROC = PAR.RPROC;
end



jobinfo = [] ;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load  trivial_regions, blocks, splits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trivial_regions = init_regions(fn_genome_config);


if ~fexist(fn_training_blocks)
  error('training block file does not exist; create it first with prepare_regions')
else
  blocks = load_struct(fn_training_blocks, 'blocks') ;
  load(fn_training_blocks, 'split') ;
  assert(length(fieldnames(split))==num_splits);
end


generate_examples=1;
if fexist([fn_candsites 'PAR.mat'])
  fprintf(1,'files with candidate pos already exists\n')
  Old = load([fn_candsites 'PAR.mat'],'P');
  diff_fields = compare_PAR(Content,Old.P.Content) ;
  if ismember('Content',Content_name)
    warning('label file exists, Content parameters are differnt!')  
    % fprintf('label file exists, Content parameters are differnt!')  
  end
  %[temp, num_files] = unix(sprintf('ls %scontig_*.pos | wc -l',fn_candsites));
  d=dir(sprintf('%scontig_*.pos',fn_candsites));
  num_files = length(d) ;

  if (num_files==length(trivial_regions))
    generate_examples=0;
    fprintf(1,'all cand files found; no label generation\n')
  else
    fprintf(1,'%i/%i cand files missing; start label generation\n',length(trivial_regions)-num_files,length(trivial_regions))
  end
  clear Old
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generating .pos .label files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if generate_examples
  fprintf(1,'generating training examples\n')
  jobinfo = rproc_empty(0) ;
  
  P.fn_genes_contents = fn_genes_contents ;
  P.fn_candsites = fn_candsites;
  P.label_source = label_source;
  P.Content = Content; 
  P.info = info;
  
  D.P = P ;
  for j=1:length(trivial_regions)
    D.t_region = trivial_regions(j);
    %keyboard
    filename = sprintf('%scontig_%i%s',fn_candsites,trivial_regions(j).chr_num,trivial_regions(j).strand);
    if fexist(filename)
      continue
    end
    if run_locally
      save_label_files_starter(D);
    else
      [mem_req, time_req, opts] = rproc_memtime_policy('save_label_files_starter',  0, RPROC.options) ;%RPROC.MEMREQ ;

      jobinfo(j) = rproc('save_label_files_starter', D, mem_req, opts, time_req);
      % save_label_files_starter(D)
    end;
    % save_label_files(t_region,PAR,genes)
  end
  if ~isempty(jobinfo),
    [jobinfo_1, num_crashed]=rproc_wait(jobinfo, 10, 1, -1);
  else
    num_crashed = 0 ;
    jobinfo_1 = rproc_empty(0) ;
  end ;
  %[temp, num_files] = unix(sprintf('ls %scontig_*.pos | wc -l',fn_candsites));
  d=dir(sprintf('%scontig_*.pos',fn_candsites));
  num_files = length(d) ;
  if num_files~=length(trivial_regions)
    error('number of files (%i) not equal to number of trivial_regions (%i)',num_files, length(trivial_regions));
  end
  if num_crashed>0, 
    num_crashed
    %return
  end 
  if ~isempty(jobinfo_1)
    rproc_cleanup(jobinfo_1) ;
  end ;
end
clear P

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generating example files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modes = {'train','test'};
for mode=modes
  [filename, exists]=gen_example_filename(Content,fn_filter_settings,fn_examples,mode{1},num_splits);
  if ~exists
        
    P.Content = Content;
    P.num_splits = num_splits;
    P.fn_training_blocks = fn_training_blocks;
    P.fn_candsites = fn_candsites;
    P.fn_filter_settings = fn_filter_settings ;
    P.fn_examples = fn_examples;
    P.fn_example_statistics = fn_example_statistics;

    if ~run_locally
      D.blocks = blocks;
      D.split = split;
      D.P = P;

      [mem_req, time_req, opts] = rproc_memtime_policy('save_content_examples',  0, RPROC.options) ;%RPROC.MEMREQ ;
      
      jobinfo = rproc('save_content_examples', D, mem_req, opts, time_req);

      [jobinfo_1, num_crashed]=rproc_wait(jobinfo, 60, 1, -1);
    else
      save_content_examples(blocks,split,P) ; num_crashed = 0 ;
    end
    if num_crashed>0
      return
    end
  end
end

