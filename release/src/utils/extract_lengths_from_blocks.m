function lengths=extract_lengths_from_blocks(blocks, type, field)

if nargin==0
	fprintf('usage: lengths=extract_lengths_from_blocks(blocks, type [, field]);\n')
end

if nargin<3
	field = 'truth';
end

lengths = zeros(1, length(blocks)*10);
lencount = 0;

for j = 1:length(blocks)	
	segments = blocks(j).(field).segments;
	starts = segments(segments(:, 3)==type, 1);
	stops = segments(segments(:, 3)==type, 2);
	len = stops - starts +1;
	if lencount+length(len)>length(lengths)
		lengths = [lengths zeros(1, length(lengths))];
	end
	lengths(lencount+1:lencount+length(len)) = len;
	lencount = lencount+length(len);
end

lengths(lencount+1:end) = [];

