function genes = load_genes(filename)

if exist(filename, 'dir')==7
  fn_anno = sprintf('%s/genes.mat',filename);
  %% for octave the whos '-file' option does not exist
  try
	load(fn_anno, 'genes')
	return
  end
elseif exist(filename, 'file')==2
  fn_anno = filename;
else
	error('file not found');
end

if ~isempty(strfind(license, 'GNU'))
	% running octave
	l = load(fn_anno, 'genes');
	genes = l.genes;
	return
end
xx = whos('-file', fn_anno);

for j = 1:length(xx)
	if strcmp(xx(j).name,'genes')
		load(fn_anno, 'genes')
		return
	elseif strcmp(xx(j).name,'all_genes')
		load(fn_anno, 'all_genes')
		genes = all_genes;
		return
	end
end

fprintf('filename: %s\n', filename);
error('variable not found');
