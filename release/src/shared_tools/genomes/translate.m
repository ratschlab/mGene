function aseq = trans(nseq)

aseq = zeros(1,floor(length(nseq)/3)) ;
j=0 ;
for i=1:3:length(nseq)-3+1
  j=j+1 ;
  aseq(j) = trans3(nseq(i:i+2)) ;
end ;
aseq = char(aseq) ;


function a=trans3(nuc)

if ((nuc(1)=='n') || (nuc(2)=='n') || (nuc(3)=='n')) 
  %a='?'; % this is not very common
  a='X'; 
  return ;
end ;

if (strncmp(nuc,'ttt',3)==1), a='F' ;
elseif (strncmp(nuc,'ttc',3)==1) a='F' ;
elseif (strncmp(nuc,'tta',3)==1) a='L' ;
elseif (strncmp(nuc,'ttg',3)==1) a='L' ;
elseif (strncmp(nuc,'ttg',3)==1) a='L' ;
elseif (strncmp(nuc,'ctt',3)==1) a='L' ;
elseif (strncmp(nuc,'ctc',3)==1) a='L' ;
elseif (strncmp(nuc,'cta',3)==1) a='L' ;
elseif (strncmp(nuc,'ctg',3)==1) a='L' ;
elseif (strncmp(nuc,'att',3)==1) a='I' ;
elseif (strncmp(nuc,'atc',3)==1) a='I' ;
elseif (strncmp(nuc,'ata',3)==1) a='I' ;
elseif (strncmp(nuc,'atg',3)==1) a='M' ;
elseif (strncmp(nuc,'gtt',3)==1) a='V' ;
elseif (strncmp(nuc,'gtc',3)==1) a='V' ;
elseif (strncmp(nuc,'gta',3)==1) a='V' ;
elseif (strncmp(nuc,'gtg',3)==1) a='V' ;
elseif (strncmp(nuc,'tct',3)==1) a='S' ;
elseif (strncmp(nuc,'tcc',3)==1) a='S' ;
elseif (strncmp(nuc,'tca',3)==1) a='S' ;
elseif (strncmp(nuc,'tcg',3)==1) a='S' ;
elseif (strncmp(nuc,'cct',3)==1) a='P' ;
elseif (strncmp(nuc,'ccc',3)==1) a='P' ;
elseif (strncmp(nuc,'cca',3)==1) a='P' ;
elseif (strncmp(nuc,'ccg',3)==1) a='P' ;
elseif (strncmp(nuc,'ac',2)==1) a='T' ;
elseif (strncmp(nuc,'gc',2)==1) a='A' ;
elseif (strncmp(nuc,'tat',3)==1) a='Y' ;
elseif (strncmp(nuc,'tac',3)==1) a='Y' ;
elseif (strncmp(nuc,'taa',3)==1) a='*' ;
elseif (strncmp(nuc,'tag',3)==1) a='*' ;
elseif (strncmp(nuc,'tga',3)==1) a='*' ;
elseif (strncmp(nuc,'cat',3)==1) a='H' ;
elseif (strncmp(nuc,'cac',3)==1) a='H' ;
elseif (strncmp(nuc,'caa',3)==1) a='Q' ;
elseif (strncmp(nuc,'cag',3)==1) a='Q' ;
elseif (strncmp(nuc,'aat',3)==1) a='N' ;
elseif (strncmp(nuc,'aac',3)==1) a='N' ;
elseif (strncmp(nuc,'aaa',3)==1) a='K' ;
elseif (strncmp(nuc,'aag',3)==1) a='K' ;
elseif (strncmp(nuc,'gat',3)==1) a='D' ;
elseif (strncmp(nuc,'gac',3)==1) a='D' ;
elseif (strncmp(nuc,'tgt',3)==1) a='C' ;
elseif (strncmp(nuc,'tgc',3)==1) a='C' ;
elseif (strncmp(nuc,'tgg',3)==1) a='W' ;
elseif (strncmp(nuc,'cg',2)==1) a='R' ;
elseif (strncmp(nuc,'aga',3)==1) a='R' ;
elseif (strncmp(nuc,'agg',3)==1) a='R' ;
elseif (strncmp(nuc,'agt',3)==1) a='S' ;
elseif (strncmp(nuc,'agc',3)==1) a='S' ;
elseif (strncmp(nuc,'gg',2)==1) a='G' ;
elseif (strncmp(nuc,'gaa',3)==1) a='E' ;
elseif (strncmp(nuc,'gag',3)==1) a='E' ;
else
  fprintf('unknown %c%c%c\n', nuc(1), nuc(2), nuc(3)) ;
  error(' ???') ;
end ;

