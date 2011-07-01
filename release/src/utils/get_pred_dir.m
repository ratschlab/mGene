function pred_dir = get_pred_dir(organism, experiment, cleave)


if 0%~isempty(strmatch(organism, {'c_brenneri', 'c_briggsae', 'c_elegans', 'c_japonica', 'c_remanei', 'p_pacificus'}, 'exact')),
	base_dir = '/fml/ag-raetsch/nobackup/projects/five_nematodes';	
else
	base_dir = '/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions';
end


if cleave==0
	pred_dir = sprintf('%s/%s/lsl/%s/output/genome_wide_predictions/', base_dir, organism, experiment);
else
	pred_dir = sprintf('%s/%s/lsl/%s/output_cleave_%i/genome_wide_predictions/', base_dir, organism, experiment, cleave);
end

organism
[a b] = unix(sprintf('ls %s', organism))

if exist(organism, 'dir')==7
	% this is a hack to make this run for the mGeneToolbox
	pred_dir = fullfile(organism, 'genome_wide_predictions');
end
