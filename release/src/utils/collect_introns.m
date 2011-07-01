function all_introns = collect_introns(genes, coding, transcripts)

if nargin==1
	coding = 0;
end

all_introns = zeros(length(genes)*10, 2); 

cnt = 0;
for j = 1:length(genes)
	if nargin<3
		transcripts = 1:length(genes(j).exons);
	end	
	for k = transcripts
		if coding
			exons = genes(j).cds_exons{k};
		else
			exons = genes(j).exons{k};
		end
		if size(exons, 1)>1
			introns = [exons(1:end-1, 2) exons(2:end, 1)];
			num = size(introns, 1);
			all_introns(cnt+1:cnt+num, :) = introns;
			cnt = cnt+num; 
		end
	end
end
all_introns(cnt+1:end, :) = [];
