function s = flip_strand(s)

if s=='+'
	s = '-';
elseif s=='-'
	s = '+';
end
