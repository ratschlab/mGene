function create_directories(FN,organism,verb,signal)

% create_directories(PAR)
%   
% function creates all necessary directories for a gene finding
% experiment as defined by PAR.
  
  
if nargin<3
  verb=0;
end

if  nargin>3
  signalnames{1} = signal;
else
  signalnames = fieldnames(FN.input_sig);
end


%%%%%%%%%%%%%%%
%%% GENERAL INPUT FILES 
%%%%%%%%%%%%%%%

exp_base_in = sprintf('%s/%s/%s',FN.input_dir,organism.name,FN.dir_name);

input_dir = sprintf('%s/orig_data/',exp_base_in);
if verb
  fprintf(1, 'dir for original data: \n%s\n\n', input_dir)
end
nice_mkdir(input_dir);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTPUT DIRs 
%%%%%%%%%%%%%%%%%%%%%%%%%

exp_base_out = sprintf('%s/%s/%s',FN.output_dir,organism.name,FN.dir_name);
output_dir = sprintf('%s/genome_data/',exp_base_out);
if verb
  fprintf(1, 'dir for genes, split and regions: \n%s\n\n', output_dir)
end
nice_mkdir(output_dir);

nice_mkdir(FN.exp_dir);


%%%%%%%%%%%%%%%
%%% SIGNAL SPECIFIC Dirs
%%%%%%%%%%%%%%%


for s=1:length(signalnames)
  if ~isfield( FN.input_sig, signalnames{s})
    continue
  end
  signals = FN.input_sig.(signalnames{s});
  signal_fns = fieldnames(signals);
  for i=1:length(signal_fns)
    name = signals.(signal_fns{i}) ;
    iii = max(strfind(name,'/'));
    if ~isempty(iii)
      nice_mkdir(name(1:iii));
    end 
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIGNAL SPECIFIC OUTPUT DIRs 
%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4
  signalnames = fieldnames(FN.output_sig);
end
  
for s=1:length(signalnames)
  if ~isfield( FN.output_sig, signalnames{s})
    continue
  end
  signals = FN.output_sig.(signalnames{s});
  signal_fns = fieldnames(signals);
  for i=1:length(signal_fns)
    name = signals.(signal_fns{i}) ;
    if isempty(strmatch(FN.output_dir,name)) &verb
      fprintf(['incomplete signal output file name, signal output directories not generated; '...
               'run config for %s\n'],signalnames{s})
      break
    end
    iii = max(strfind(name,'/'));
    if ~isempty(iii)
      nice_mkdir(name(1:iii));
    end 
  end
end

%%%%%%%%%%%%%%%
%%% CONTENT SPECIFIC Dirs
%%%%%%%%%%%%%%%

if isfield(FN, 'input_cont')
  if nargin<4
    signalnames = fieldnames(FN.input_cont);
  end
  for s=1:length(signalnames)
    if ~isfield( FN.input_cont, signalnames{s})
      continue
    end
    signals = FN.input_cont.(signalnames{s});
    signal_fns = fieldnames(signals);
    for i=1:length(signal_fns)
      name = signals.(signal_fns{i}) ;
      iii = max(strfind(name,'/'));
      if ~isempty(iii)
        nice_mkdir(name(1:iii));
      end 
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONTENT SPECIFIC OUTPUT DIRs 
%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(FN, 'output_cont')
  if nargin<4
    signalnames = fieldnames(FN.output_cont);
  end
  for s=1:length(signalnames)
    if ~isfield( FN.output_cont, signalnames{s})
      continue
    end
    signals = FN.output_cont.(signalnames{s});
    signal_fns = fieldnames(signals);
    for i=1:length(signal_fns)
      name = signals.(signal_fns{i}) ;
      if isempty(strmatch(FN.output_dir,name)) &verb
        fprintf(['incomplete signal output file name, signal output directories not generated; '...
                 'run create_filenames first for\n'])
        break
      end
      iii = max(strfind(name,'/'));
      if ~isempty(iii)
        nice_mkdir(name(1:iii));
      end 
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%% LSL SPECIFIC
%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4
  
  if isfield(FN, 'output_lsl')
    fns = fieldnames(FN.output_lsl);
    for i=1:length(fns)
      name = FN.output_lsl.(fns{i}) ;
      if isempty(strmatch(FN.output_dir,name)) &verb
        fprintf('incomplete signal output file name, signal output directories not generated;')
        break
      end
      iii = max(strfind(name,'/'));
      if ~isempty(iii)
        nice_mkdir(name(1:iii));
      end 
    end
  end
  if isfield(FN, 'input_lsl')
    fns = fieldnames(FN.input_lsl);
    for i=1:length(fns)
      if ~iscell(FN.input_lsl.(fns{i}))
        name = FN.input_lsl.(fns{i}) ;
        if isempty(strmatch(FN.input_dir,name)) &verb
          fprintf('incomplete signal input file name, signal input directories not generated;')
          break
        end
        iii = max(strfind(name,'/'));
        if ~isempty(iii)
          nice_mkdir(name(1:iii));
        end 
      else
        for j=1:length(FN.input_lsl.(fns{i}))
          name = FN.input_lsl.(fns{i}){j} ;
          if isempty(strmatch(FN.input_dir,name)) &verb
            fprintf('incomplete signal input file name, signal input directories not generated;')
            break
          end
          iii = max(strfind(name,'/'));
          if ~isempty(iii)
            nice_mkdir(name(1:iii));
          end 
        end
      end
    end
  end
  
  %fn = sprintf('%s/%s/current_predictions', FN.input_dir, organism.name);
  %nice_mkdir(fn);

end
