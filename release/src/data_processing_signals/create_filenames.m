function FN = create_filenames(FN,organism,Signals,Contents,verb)

%
% FN requires fields :
% 	- input_dirC_
% 	- output_dir
%	- exp_name
%	- FN.fn_genes_train
%	- date
%	- user
%	- organism
%

if nargin<2
  verb=0;
end

if ~isfield(FN,'dir_name_ext')
  dir_name_ext = '';
else
  dir_name_ext = [ '_', FN.dir_name_ext ];
end;

if isempty(FN.dir_name)
  FN.dir_name ='';
end
if isempty(FN.exp_name)
  FN.exp_name = '';
end
if isempty(FN.lsl_exp_name)
  FN.lsl_exp_name = '';
end

assert(isfield(FN,'exp_name'))
assert(isfield(FN,'dir_name'))

%%%%%%%%%%%%%%%
%%% GENERAL INPUT GENOME FILES 
%%%%%%%%%%%%%%%
exp_base_in = sprintf('%s/%s/%s/%s',FN.genome_dir,organism.name,organism.release,organism.genebuilt);

input_fns = fieldnames(FN.genome);
for i=1:length(input_fns)
  old_name = getfield(FN.genome,input_fns{i});
  if iscell(old_name)
    new_name ={};
    for j=1:length(old_name)
      if isempty(strmatch(FN.genome_dir,old_name{j}))
        new_name{j}  = sprintf('%s/%s',exp_base_in,old_name{j} );
      else
        new_name{j}  = old_name{j} ;
      end
    end
  elseif length(strfind(old_name,'/'))==1
    exp_base_in_ = sprintf('%s/%s/%s',FN.genome_dir,organism.name,organism.release);
    new_name  = sprintf('%s/%s',exp_base_in_,old_name );
  elseif isempty(strmatch(FN.genome_dir,old_name)) 
    new_name  = sprintf('%s/%s',exp_base_in,old_name );
  else
    new_name  = old_name ;
  end
  FN.genome = setfield(FN.genome,input_fns{i},new_name);
end  
if verb
  fprintf(1, 'genebuild from dir: \n%s\n\n', exp_base_in)
end

%%%%%%%%%%%%%%%
%%% GENERAL INPUT FILES 
%%%%%%%%%%%%%%%

exp_base_in = sprintf('%s/%s/%s',FN.input_dir,organism.name,FN.dir_name);

input_fns = fieldnames(FN.input);
for i=1:length(input_fns)
  old_name = getfield(FN.input,input_fns{i});
  if iscell(old_name)
    new_name ={};
    for j=1:length(old_name)
      if isempty(strmatch(FN.input_dir,old_name{j}))
        new_name{j}  = sprintf('%s/orig_data/%s',exp_base_in,old_name{j} );
      else
        new_name{j}  = old_name{j} ;
      end
    end
  else
    if isempty(strmatch(FN.input_dir,old_name)) 
      new_name  = sprintf('%s/orig_data/%s',exp_base_in,old_name );
    else
      new_name  = old_name ;
    end
  end
  FN.input = setfield(FN.input,input_fns{i},new_name);
end  
input_dir=sprintf('%s/orig_data/',exp_base_in);
if verb
  fprintf(1, 'dir for original data: \n%s\n\n', input_dir)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTPUT DIRs AND FNs  
%%%%%%%%%%%%%%%%%%%%%%%%%
exp_base_out = sprintf('%s/%s/%s',FN.output_dir,organism.name,FN.dir_name);

output_fns = fieldnames(FN.output);
for i=1:length(output_fns)
  old_name = getfield(FN.output,output_fns{i});
  if isempty(old_name)
    old_name = '';
  end 
  if iscell(old_name)
    new_name ={};
    for j=1:length(old_name)
      if isempty(strmatch(FN.output_dir,old_name{j})) 
        new_name{j}  = sprintf('%s/genome_data/%s',exp_base_out,old_name{j} );
      else
        new_name{j}  = old_name{j} ;
      end
    end
  else
    if  isempty(strmatch(FN.output_dir,old_name)) 
      new_name  = sprintf('%s/genome_data/%s',exp_base_out,old_name );
    else
      new_name  = old_name ;
    end
  end
  FN.output = setfield(FN.output,output_fns{i},new_name);
end  

output_dir=sprintf('%s/genome_data/',exp_base_out);
if verb
  fprintf(1, 'dir for genes, split and regions: \n%s\n\n', output_dir)
end
if ~isfield( FN,'exp_dir')||isempty(FN.exp_dir)|| isempty(strmatch(FN.input_dir,FN.exp_dir)) 
  FN.exp_dir = sprintf('%s/%s',exp_base_out,FN.exp_name);
end

%%%%%%%%%%%%%%%
%%% SIGNAL SPECIFIC INPUT FILES 
%%%%%%%%%%%%%%%

signalnames = fieldnames(FN.input_sig);
for s=1:length(signalnames)
  signals = getfield(FN.input_sig,signalnames{s});
  signal_fns = fieldnames(signals);
  for i=1:length(signal_fns)
    old_name = getfield(signals,signal_fns{i}) ;
    if isempty(old_name)
       old_name = '';
    end 
    if isempty(strmatch(FN.input_dir,old_name)) 
      new_name  = sprintf('%s/%s/%s/%s',exp_base_in,FN.exp_name,dir_name_ext,old_name );
      % new_name  = sprintf('%s/%s/sensors/%s%s/%s',exp_base_in,FN.exp_name,signalnames{s},dir_name_ext,old_name );
    else
      new_name  = old_name ; 
    end
    signals = setfield(signals,signal_fns{i},new_name);
  end 
  FN.input_sig = setfield(FN.input_sig,signalnames{s},signals);
  
  % sig_input = sprintf('%s/%s/sensors/%s%s',exp_base_in,FN.exp_name,signalnames{s},dir_name_ext,old_name);
  if verb
    fprintf(1, 'dir for label_files and examples: \n%s\n\n', new_name)
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIGNAL SPECIFIC OUTPUT DIRs AND FNs  
%%%%%%%%%%%%%%%%%%%%%%%%%



%output_fns = fieldnames(FN.output);
%for i=1:length(output_fns)
%  old_name = getfield(FN.output,output_fns{i});
%  new_name  = sprintf('%s/%s_%s/%s',FN.output_dir,FN.exp_name,...
%                      organism.name,old_name );
%  FN.output = setfield(FN.output,output_fns{i},new_name);
%end  

signalnames = fieldnames(FN.output_sig);
for s=1:length(signalnames)
  signals = FN.output_sig.(signalnames{s});
  signal_fns = fieldnames(signals);
  method = Signals.(signalnames{s}).method;
  for i=1:length(signal_fns)
    old_name = getfield(signals,signal_fns{i}) ;
    dir_name = sprintf('%s/%s/%s',exp_base_out,FN.exp_name,dir_name_ext);
    % dir_name = sprintf('%s/%s/sensors/%s%s',exp_base_out,FN.exp_name,signalnames{s},dir_name_ext);
    if isfield(FN, 'custom_output')
      dir_name = sprintf('%s/%s',dir_name,FN.custom_output);
    else
      dir_name = sprintf('%s/output_%s',dir_name,method.name);
      if isequal(method.name, 'SVM')
        if length(method.kernels)==1&&~isempty(method.kernels{1})
          if isfield(method.kernels{1},'name')
            kernel_names = method.kernels{1}.name; 
          else
            kernel_names = method.kernels{1}.kernel_name;
          end
        else
          kernel_names = 'COMBINED';
        end
        dir_name = sprintf('%s%s',dir_name,kernel_names);
      end
    end
    if isempty(strmatch(FN.output_dir,old_name)) 
      new_name  = sprintf('%s/%s',dir_name, old_name );
    else
      new_name = old_name;
    end
    signals = setfield(signals,signal_fns{i},new_name);
  end
  FN.output_sig = setfield( FN.output_sig, signalnames{s}, signals );
  if verb
    fprintf(1, 'dir for training and test data: \n%s\n\n', dir_name) 
  end
end


%%%%%%%%%%%%%%%
%%% CONTENT SPECIFIC INPUT FILES 
%%%%%%%%%%%%%%%

if isfield(FN, 'input_cont')
  cont_names = fieldnames(FN.input_cont);
  for s=1:length(cont_names)
    content = getfield(FN.input_cont,cont_names{s});
    cont_fns = fieldnames(content);
    for i=1:length(cont_fns)
      old_name = getfield(content,cont_fns{i}) ;
      if isempty(strmatch(FN.input_dir,old_name)) 
        new_name  = sprintf('%s/%s/%s/%s',exp_base_in,FN.exp_name,dir_name_ext,old_name );
        % new_name  = sprintf('%s/%s/contents/%s%s/%s',exp_base_in,FN.exp_name, cont_names{s},dir_name_ext,old_name );
      else
        new_name  = old_name ; 
      end
      content = setfield(content,cont_fns{i},new_name);
    end 
    FN.input_cont = setfield(FN.input_cont,cont_names{s},content);
  end

  %cont_input = sprintf('%s/%s/contents/%s',exp_base_in,FN.exp_name,signal_dir_name);
  if verb
    fprintf(1, 'dir for content specific input: \n%s\n\n', cont_input)
  end
end

%%%%%%%%%%%%%%%
%%% CONTENT SPECIFIC OUTPUT FILES 
%%%%%%%%%%%%%%%

if isfield(FN, 'output_cont')
  cont_names = fieldnames(FN.output_cont);
  for s=1:length(cont_names)
    content = getfield(FN.output_cont,cont_names{s});
    cont_fns = fieldnames(content);
    method = Contents.(cont_names{s}).method;
    for i=1:length(cont_fns)
      old_name = getfield(content,cont_fns{i}) ;
      dir_name = sprintf('%s/%s/%s',exp_base_out,FN.exp_name,dir_name_ext);
      if isfield(FN, 'custom_output')
        dir_name = sprintf('%s/%s',dir_name,FN.custom_output);
      else
        dir_name = sprintf('%s/output_%s',dir_name,method.name);
        if isequal(method.name, 'SVM')
          if length(method.kernels)==1&&~isempty(method.kernels{1})
            if isfield(method.kernels{1},'name')
              kernel_names = method.kernels{1}.name; 
            else
              kernel_names = method.kernels{1}.kernel_name;
            end
          else
            kernel_names = 'COMBINED';
          end
          dir_name = sprintf('%s%s',dir_name,kernel_names);
        end
      end
      if isempty(strmatch(FN.output_dir,old_name)) 
        new_name  = sprintf('%s/%s',dir_name, old_name );
      else
        new_name = old_name;
      end
      content = setfield(content,cont_fns{i},new_name);
    end 
    FN.output_cont = setfield(FN.output_cont, cont_names{s},content);
  end

  %cont_output = sprintf('%s/%s/contents/%s',exp_base_out,FN.exp_name,signal_dir_name);
  if verb
    fprintf(1, 'dir for content specific output: \n%s\n\n', cont_output)
  end
end

%%%%%%%%%%%%%%%
%%% LSL STUFF
%%%%%%%%%%%%%%%

if isfield(FN, 'input_lsl')
  fns = fieldnames(FN.input_lsl);
  FN.input_lsl.dir = sprintf('%s/%s/%s/%s/lsl/%s/', ...
			     exp_base_out,organism.name,FN.dir_name,FN.exp_name,FN.lsl_exp_name ); 
  for i=1:length(fns)
    if ~iscell(FN.input_lsl.(fns{i}))
      old_name = FN.input_lsl.(fns{i}) ;
      if isempty(strmatch(FN.input_lsl.dir,old_name))&~isequal(fns{i},'fn_old_training_blocks')...
            &~isequal(fns{i},'fn_old_test_blocks')&~isequal(fns{i},'fn_old_split')
        new_name  = sprintf('%s/%s',FN.input_lsl.dir,old_name);
      else
        new_name  = old_name ; 
      end
      FN.input_lsl.(fns{i}) = new_name;
    else
      for j=1:length(FN.input_lsl.(fns{i}))
        old_name = FN.input_lsl.(fns{i}){j} ;
        if isempty(strmatch(FN.input_lsl.dir,old_name))&~isequal(fns{i},'fn_old_training_blocks')...
              &~isequal(fns{i},'fn_old_test_blocks')&~isequal(fns{i},'fn_old_split')
          new_name  = sprintf('%s/%s',FN.input_lsl.dir,old_name);
        else
          new_name  = old_name ; 
        end
        FN.input_lsl.(fns{i}){j} = new_name;
      end    
    end
  end
end
%cont_output = sprintf('%s/%s/contents/%s',exp_base_out,FN.exp_name,signal_dir_name);
if verb
  fprintf(1, 'dir for lsl specific input: \n%s\n\n', cont_output)
end

if isfield(FN, 'output_lsl')
  fns = fieldnames(FN.output_lsl);
  FN.output_lsl.dir = sprintf('%s/%s/%s/%s/lsl/%s/', ...
			      exp_base_out, organism.name,FN.dir_name,FN.exp_name,FN.lsl_exp_name ); 
  for i=1:length(fns)
    old_name = FN.output_lsl.(fns{i}) ;
    if isempty(strmatch(FN.output_lsl.dir,old_name))
      new_name  = sprintf('%s/%s',FN.output_lsl.dir,old_name);
    else
      new_name  = old_name ; 
    end
    FN.output_lsl.(fns{i}) = new_name;
  end 
end
%cont_output = sprintf('%s/%s/contents/%s',exp_base_out,FN.exp_name,signal_dir_name);
if verb
  fprintf(1, 'dir for lsl specific output: \n%s\n\n', cont_output)
end

