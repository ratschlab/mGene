function ok = check_predictions_helper(P)
% ok = check_predictions(P)

fn_genome_config = P.fn_genome_config;
chrom = P.chrom;
strand = P.strand;
fn_pred = P.fn_pred;
resolution = P.resolution;
Signal_name = P.Signal_name;

genome_info = init_genome(fn_genome_config);
S = dir(genome_info.flat_fnames{chrom}) ;
contig_length = S.bytes;
CHR_NAME = genome_info.contig_names{chrom};
fprintf('contig: %i%s\n', chrom,strand);

%% GET ALL PREDICTIONS from mat file
filename = sprintf('%scontig_%i%s_all.mat', fn_pred, ...
                   chrom, strand);
load(filename);
ok=1;


% always save global coordinates on forward strand
if isequal(strand, '-'),
  pos = contig_length - pos + 1;
end
if ~(all(pos>=1));
  fprintf('global positions wrong: pos<1')
  ok=-1;
end

[pos sort_idx] = sort(pos);
output = output(sort_idx);

if exist('Conf', 'var') &&exist('Conf_cum', 'var')
  Conf = Conf(sort_idx);
  Conf_cum = Conf_cum(sort_idx);

  if ~(all(Conf>=0-1e-8&Conf<=1+1e-8));
    fprintf('Confs out of range')
  end
end

trivial_regions = init_regions(fn_genome_config);

idx  = find([trivial_regions.strand] == strand&[trivial_regions.chr_num] == chrom);
reg = trivial_regions(idx);

PAR.FN.input.fn_genome_config = fn_genome_config ;
PAR.FN.output_sig.(Signal_name).fn_pred = fn_pred ;

if exist('Conf', 'var') &&exist('Conf_cum','var')
  fields = {'output','Conf','Conf_cum'};
else
  fields = {'output'} ;
end
reg = add_signals2blocks(reg, PAR, {Signal_name}, '', 0,fields);
 
[tmp,idx1,idx2] = intersect(reg.Signals.(Signal_name).contig_pos,pos);
num_pos = max(length(pos),length(reg.Signals.(Signal_name)));
if length(tmp)<num_pos
  fprintf('WARNING: found only %i of %i  postions (%i %%)\n',length(tmp),num_pos,round((length(tmp))/num_pos*100));
  fprintf('specified resolution:%i\n',resolution);
end

if isfield(reg.Signals.(Signal_name),'Conf_cum') && ~(all(abs(reg.Signals.(Signal_name).Conf_cum(idx1)-Conf_cum(idx2)')<1e-3))
  fprintf('something wrong with file')
  ok=-1;
end
