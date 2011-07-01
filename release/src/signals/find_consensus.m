function POS=  find_consensus(seq,Signal,strand)
%  POS=  find_consensus(seq,Signal,strand) 

consensus = Signal.consensus ; 
lwin_big = Signal.lwin_big;
rwin_big = Signal.rwin_big;
offset = Signal.offset;



POS = zeros(1,floor(length(seq)/10));
num = 0;
seq=upper(seq);
for m = 1:length(consensus)  
  cons = upper(consensus{m});
  pos_tmp = findstr(seq(lwin_big+1:end-rwin_big),cons)+lwin_big;
  POS(num+[1:length(pos_tmp)]) = pos_tmp;
  num = num + length(pos_tmp);
end
POS(num+1:end)=[];
clear num 
POS = POS + offset;
%if strand=='-' & strcmp(Signal.name,'polya')
%  POS = -POS + length(seq)+1-length(cons)+1;
%elseif strand=='-' 
if strand=='-'
  POS = -POS + length(seq)+1;
end
POS = unique(sort(POS));
