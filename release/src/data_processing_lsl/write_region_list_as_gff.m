function write_region_list_as_gff(blocks, file, append_flag, info)

if nargin==4&&append_flag
	[fd msg] = fopen(file, 'a+');%append
else
	[fd msg] = fopen(file, 'w+');%discard
end

if fd<1, error('could not open file %s: %s', file, msg); end

for j=1:length(blocks)
	chr = blocks(j).chr;
	if isfield(blocks, 'start')
		start = blocks(j).start;
		stop = blocks(j).stop;
	else
		start = blocks(j).block_start;
		stop = blocks(j).block_end;
	end
	strand = blocks(j).strand;
	if nargin==4
		fprintf(fd,'%s\t%i\t%i\t%s\t%s\n',chr,start,stop,strand,info);
	else
		fprintf(fd,'%s\t%i\t%i\t%s\n',chr,start,stop,strand);
	end
end

fclose(fd);
