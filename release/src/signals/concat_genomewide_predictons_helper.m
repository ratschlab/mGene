function [dummy1, dummy2]=concat_genomewide_predictons_helper(P)
% concat_genomewide_predictons_helper(P)

try
  pause('on')
  warning('off', 'MATLAB:typeaheadBufferOverflow');
catch
  % works only for matlab
end


dummy1=[] ;
dummy2=[] ;

fn_pred_all = P.fn_pred_all ;
fn_pred_split = P.fn_pred_split;
num_splits = P.num_splits;
avg_all_svms = P.avg_all_svms;
subtract_diff = P.subtract_diff;

% for H_sapiens_nohomologs:
% prefix = '/fml/ag-raetsch/share/projects/genefinding/H_sapiens_nohomologs/sensors/don/output/pred';
% fn_pred_all = sprintf('%scontig_%i%s_all.mat',prefix,c,s);

contig_pos = [];
contig_out = [];
contig_conf= [];
contig_cnfc= [];
contig_svm = [];
for p=1:num_splits
  fn_pred = sprintf('%s%i.mat', fn_pred_split, p);
  % for H_sapiens_nohomologs:
  % prefix = '/fml/ag-raetsch/share/projects/genefinding/H_sapiens_nohomologs/sensors/don/output/pred';
  % fn_pred = sprintf('%scontig_%i%s%s_%i.mat',prefix,c,s,PAR.Signal_name,p);

  for kk=1:10
    try
      clear pos output Conf Conf_cum
      load(fn_pred, 'inventory');
      %inventory
      for j=1:length(inventory)
         if kk>1,
           fprintf('trying to load variable "%s" from %s (trial %i)\n', inventory{j}, fn_pred, kk);
         end ;
	 load(fn_pred,inventory{j})
	 assert(exist(inventory{j},'var')==1);	 
      end
      break ;
    catch
      if kk==10      
        fprintf('giving up waiting for file %s\n', fn_pred)
      else
        fprintf('waiting for file: %s \n', fn_pred);    
        pause(30)
      end	
    end
  end

  if ~avg_all_svms
    if ~isempty(intersect(contig_pos, pos))
      warning('positions found in multiple splits');
      %%hack
      [tmp,idx1,idx2] = intersect(contig_pos, pos);
      pos(idx2)=[];
      output(idx2)=[];
      if exist('Conf', 'var') && exist('Conf_cum','var')
        Conf(idx2)=[];
        Conf_cum(idx2)=[];
      else
        warning('no confs found')
      end
      %%hack
    end
    % assert(isempty(intersect(contig_pos, pos)));
    contig_pos = [contig_pos pos];
    contig_out = [contig_out output];
    contig_svm = [contig_svm repmat(p,size(pos))];
    if exist('Conf', 'var') && exist('Conf_cum','var')
      contig_conf = [contig_conf Conf];
      contig_cnfc = [contig_cnfc Conf_cum];
    end
    
  else
    [tmp,idx1,idx2] = intersect(contig_pos, pos);
    assert(length(idx1)==length(contig_pos))
    if isempty(contig_pos)
      contig_pos = pos;
      contig_out = output;
      if exist('Conf', 'var') && exist('Conf_cum', 'var')
        contig_conf = Conf(idx2);
        contig_cnfc = Conf_cum(idx2) ;
      end
    else
      contig_out(idx1) = contig_out(idx1)+output(idx2);
      if exist('Conf', 'var') && exist('Conf_cum','var')
        contig_conf(idx1) = contig_conf(idx1) + Conf(idx2) ;
        contig_cnfc(idx1) = contig_cnfc(idx1) + Conf_cum(idx2) ;
      end
    end
  end
  clear output pos Conf Conf_cum
end 
[pos sort_idx] = sort(contig_pos);
output = contig_out(sort_idx);
svm = contig_svm(sort_idx);

if avg_all_svms
  output = output./num_splits;
end

% for contents a bit of a hack
if subtract_diff
  ii = find(svm(1:end-1)~=svm(2:end));
  fprintf('number of blocks: %i\n',length(ii));
  for j=ii
    output(j+1:end) = output(j+1:end) - (output(j+1)-output(j));
  end
end

if isequal(size(contig_conf), size(contig_pos)) && isequal(size(contig_cnfc), size(contig_pos))
    
  Conf = contig_conf(sort_idx);
  Conf_cum = contig_cnfc(sort_idx);
  if avg_all_svms
    Conf = Conf./num_splits;
    Conf_cum = Conf_cum./num_splits;
  end
  assert(all(Conf>=0-1e-8&Conf<=1+1e-8));
  assert(all(Conf_cum>=0-1e-8&Conf_cum<=1+1e-8));
  fprintf('saving output, pos, Conf and Conf_cum\n')
  fprintf('save: %s\n',fn_pred_all)
  inventory = {'output', 'pos', 'Conf', 'Conf_cum'};
  save(fn_pred_all, '-V7','inventory', 'output', 'pos', 'Conf', 'Conf_cum')
else
  assert(isempty(contig_conf));
  assert(isempty(contig_cnfc));
  inventory = {'output', 'pos'};
  fprintf('saving output and pos without Conf and Conf_cum\n')
  fprintf('save: %s\n',fn_pred_all)
  save(fn_pred_all, '-V7','inventory', 'output', 'pos')
end
