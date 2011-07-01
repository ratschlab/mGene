function blocks = parse_region_list(gff_file);

[fd msg] = fopen(gff_file, 'r');

if fd<1, error('could not open file %s: %s',gff_file, msg); end

blocks 			= struct;
blocks.strand 		= [];
blocks.chr 		= [];
blocks.chr_num 		= [];
blocks.start 		= [];
blocks.stop 		= [];
blocks.gene_ids 	= [];
blocks.paralogs 	= [];
blocks.is_alt 		= [];
blocks.gene_complete 	= [];
blocks.config 		= [];
blocks.split 		= [];


count = 0;
while ~feof(fd)
	Line = fgetl(fd);
	if Line(1)=='#',
		continue ;
	end
	count = count+1;
	elems = separate(Line);
	if isempty(elems{1})
		error('region-list file invalid %s', gff_file);
	end
	if length(elems) < 4
		error('region-list file invalid %s: less than 4 tab delimited columns', gff_file);
	end
	start = str2num(elems{2});
	if isempty(start)
		error('region-list file invalid %s: colums 2 not a number %s', gff_file, elems{2});
	end
	stop	= str2num(elems{3});
	if isempty(stop)
		error('region-list file invalid %s: column 3 not a number %s', gff_file, elems{3});
	end


	blocks(count).strand 		= elems{4};
	blocks(count).chr 		= elems{1};
	blocks(count).chr_num 	= [];
	blocks(count).start 		= start;
	blocks(count).stop 		= stop;
	blocks(count).gene_ids	= [];
	blocks(count).paralogs	= [];
	blocks(count).is_alt 		= [];
	blocks(count).gene_complete 	= [];
	blocks(count).config 		= [];
	blocks(count).split 		= [1];

	if length(elems)>4
		if strcmp(elems{5}, 'training')
			blocks(count).train_val_split = 1;
		elseif strcmp(elems{5}, 'validation')
			blocks(count).train_val_split = -1;
		end
	end	
end

fclose(fd);
