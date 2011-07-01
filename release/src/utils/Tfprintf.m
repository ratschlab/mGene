function Tfprintf(fids, varargin)

for fid=unique(fids)
	fprintf(fid, varargin{:});
end 
