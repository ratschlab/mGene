function [dummy1, dummy2]=save_predictions_helper(P)
% save_predictions_helper(P) 
try
  pause('on')
  warning('off', 'MATLAB:typeaheadBufferOverflow');  
catch
  % works only for matlab
end


dummy1=[] ;
dummy2=[] ;

fn_genome_config = P.fn_genome_config;
chrom = P.chrom;
strand = P.strand;
fn_pred = P.fn_pred;
Signal_name = P.Signal_name;
resolution = P.resolution;
conf_cum_thresh = P.conf_cum_thresh;
do_gzip = P.do_gzip ;
save_as = P.save_as ;

genome_info = init_genome(fn_genome_config); 
S = dir(genome_info.flat_fnames{chrom}) ;
contig_length = S.bytes;
CHR_NAME = genome_info.contig_names{chrom};
fprintf('contig: %i%s\n', chrom,strand);

%----------------------------
%% GET ALL PREDICTIONS from mat file
%----------------------------

filename = sprintf('%scontig_%i%s_all.mat', fn_pred, chrom, strand) ;
for kk=1:30 % it can take easily several minutes until the file is visible
   try 
     clear pos output Conf Conf_cum
     load(filename, 'inventory');
     for j=1:length(inventory)
        if kk>1,
           fprintf('trying to load variable "%s" from %s (trial %i)\n', inventory{j}, filename, kk);
        end ;
        load(filename,inventory{j})
        assert(exist(inventory{j}, 'var')==1);	 
     end
     break ;
   catch
     if kk==30      
       fprintf('giving up waiting for file %s\n', filename)
     else
       fprintf('waiting for file: %s \n', filename);    
       pause(30)
     end	
  end
end


% always save global coordinates on forward strand
if isequal(strand, '-'),
  pos = contig_length - pos + 1;
end
assert(all(pos>=1));

[pos sort_idx] = sort(pos);
output = output(sort_idx);
if exist('Conf', 'var') && exist('Conf_cum', 'var')
  Conf = Conf(sort_idx);
  Conf_cum = Conf_cum(sort_idx);
  
  %----------------------------
  %% downsample if demanded
  %----------------------------
  
  assert(all(Conf>=0-1e-8&Conf<=1+1e-8));
  Conf(Conf<0)=0 ;
  Conf(Conf>1)=1 ;
end

side = floor(resolution/2) ;
if 0%(resolution~=1), 
  fprintf('Reducing resolution to one prediction per %i nucleotides\n', resolution) ;
  idx = 1:resolution:length(pos) ;
  cc = zeros(1,length(idx)) ;
  co = zeros(1,length(idx)) ;
  oo = zeros(1,length(idx)) ;
  pp = zeros(1,length(idx)) ;
  
  for j=1:length(idx),
    %if mod(j,1000)==0, fprintf(' %i  \r', j), end ;
    win = max(1,idx(j)-side):min(idx(j)+side,length(pos)) ;
    neighbour = find((pos(win)>=pos(idx(j))-side) & (pos(win)<pos(idx(j))+side)) ;
    %pos(win(neighbour))
    if exist('Conf', 'var') && exist('Conf_cum', 'var')
      [mv,mi] = max(Conf_cum(win(neighbour))) ;
      assert(Conf_cum(win(neighbour(mi)))==mv) ;
      cc(j) = Conf_cum(win(neighbour(mi))) ;
      co(j) = Conf(win(neighbour(mi))) ;
    else
      [mv,mi] = max(output(win(neighbour))) ;
    end
    pp(j) = pos(win(neighbour(mi))) ;
    oo(j) = output(win(neighbour(mi))) ;
  end;
  pos = pp ;
  Conf_cum = cc ;
  Conf = co ;
  output = oo ;
end ;

%----------------------------
%%check if files exist
%----------------------------
filename = sprintf('%scontig_%i%s', fn_pred, chrom, strand);
if fexist([filename '.pos']),
  fprintf('removing position file (%s). \n', filename);
  unix(sprintf('rm %s', [filename '.pos']));
end
if fexist([filename '.output']),
  fprintf('removing output file (%s)\n', filename);
  unix(sprintf('rm %s', [filename '.output']));
end
if fexist([filename '.Conf']),
  fprintf('removing Conf file (%s)\n', filename);
  unix(sprintf('rm %s', [filename '.Conf']));
end
if fexist([filename '.Conf_cum']),
  fprintf('removing Conf_cum file (%s)\n', filename);
  unix(sprintf('rm %s', [filename '.Conf_cum']));
end

%Conf=round(Conf*1000)/1000;
%output = round(output*1000)/1000;
%Conf_cum=round(Conf_cum*1000)/1000; 
if exist('Conf', 'var') && exist('Conf_cum', 'var') %&& ~all(Conf_cum==0)&&~all(Conf==0)
  save_score_pos(pos, [output; Conf; Conf_cum], filename, {'output', 'Conf', 'Conf_cum'});
else
  save_score_pos(pos, [output], filename, {'output'});
end

if isequal(Signal_name, 'transacc'),
  disp('Using acceptor predictions in combination with transacc predictions (for wiggle files only)') ;
  FN2.fn_pred = strrep(fn_pred,'transacc', 'acc') ;
  
  filename2 = sprintf('%scontig_%i%s_all.mat', FN2.fn_pred,chrom,strand);
  L=load(filename2, 'pos', 'Conf', 'Conf_cum', 'output');
  % always save global coordinates on forward strand
  if isequal(strand, '-'),
    L.pos = contig_length - L.pos + 1;
  end
  assert(all(pos>=1));
  if ~isequal(sort(L.pos),pos)
    warning('position lists not identical') ;
    setdiff(L.pos, pos)
    setdiff(pos, L.pos)
    
    [tmp,idx1,idx2]=intersect(pos,L.pos) ;
    
    L.Conf = L.Conf(idx2) ;
    L.Conf_cum = L.Conf_cum(idx2) ;
    L.output = L.output(idx2) ;
    L.pos = L.pos(idx2) ;
    
    Conf = Conf(idx1) ;
    Conf_cum = Conf_cum(idx1) ;
    output = output(idx1) ;
    pos = pos(idx1) ;
  else
    [L.pos sort_idx] = sort(L.pos);
    L.Conf = L.Conf(sort_idx);
    L.Conf_cum = L.Conf_cum(sort_idx);
    L.output = L.output(sort_idx);
  end ;
  
  Conf_cum = Conf_cum .* L.Conf_cum ;
  Conf = Conf .* L.Conf ;
  output = output + L.output ;
  clear L ;
end ;

%----------------------------
%% SAVE PREDICTIONS as SPF files
%---------------------------- 

% saving in signal prediction format (without filtering)
if save_as.spf_ascii,
  save_predictions_as_sigpred(filename, pos, output, 'output', Signal_name, CHR_NAME, strand, do_gzip);
end ;
if save_as.spf_binary,
  save_predictions_as_sigpred_bin(filename, pos, output, 'output', Signal_name, CHR_NAME, strand);
end ;

if exist('Conf', 'var') %&& ~all(Conf==0)
  if save_as.spf_ascii,
    save_predictions_as_sigpred(filename, pos, Conf, 'Conf', Signal_name, CHR_NAME, strand, do_gzip);
  end ;
  if save_as.spf_binary,
    save_predictions_as_sigpred_bin(filename, pos, Conf, 'Conf', Signal_name, CHR_NAME, strand);
  end ;
end ;

if exist('Conf_cum','var') %&&  ~all(Conf_cum==0)
  if save_as.spf_ascii,
    save_predictions_as_sigpred(filename, pos, Conf_cum, 'Conf_cum', Signal_name, CHR_NAME, strand, do_gzip);
  end ;
  if save_as.spf_binary,
    save_predictions_as_sigpred_bin(filename, pos, Conf_cum, 'Conf_cum', Signal_name, CHR_NAME, strand);
  end
end ;

%----------------------------
%% SAVE PREDICTIONS as wiggle files
%---------------------------- 

% rounding for wiggle format output
output = round(output*1000)/1000;

%----------------------------
% REMOVE PREDICTIOND BELOW THRESHOLD
%----------------------------


if exist('Conf', 'var') && exist('Conf_cum', 'var') && save_as.wiggle

  fprintf('Removing all predictions below %.2f before saving wiggle file',conf_cum_thresh) 
  idx = find(Conf > conf_cum_thresh);
  output   = output(idx);
  pos = pos(idx);
  Conf = round(Conf*1000)/1000;
  Conf_cum = round(Conf_cum*1000)/1000; 
  Conf     = Conf(idx); 
  Conf_cum = Conf_cum(idx);
  if save_as.wiggle
    save_predictions_as_wiggle(filename, pos, Conf, 'Conf', Signal_name, CHR_NAME, strand, do_gzip);
    save_predictions_as_wiggle(filename, pos, Conf_cum, 'Conf_cum', Signal_name, CHR_NAME, strand, do_gzip);
  end ;

end

if save_as.wiggle
  save_predictions_as_wiggle(filename, pos, output, 'output', Signal_name, CHR_NAME, strand, do_gzip);
end
