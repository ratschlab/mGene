function viewsplicegraph_compare(gene,exon)
% function viewsplicegraph_compare(gene,exon)
% gene is an element of the genes data structure.
% exon is a pair of coordinates

viewgraph_compare(gene.exons,gene.name, gene.transcripts, [], gene.strands, ...
	  gene.splicegraph{1},gene.splicegraph{2},exon);
return


function viewgraph_compare(truepos, name, labels,colors, strands,...
			   vertices,edges,exon)

clf;

start=inf ;
stop=-inf ;
for k=1:length(truepos)
  start=min(start,truepos{k}(1,1)-50) ;
  stop=max(stop,truepos{k}(end,2)+50) ;
end ;

cnt=0;
width=0.075 ; width2=0.03; step=0.75 ;

gridsize = round((stop-start)/5);
for i=ceil(start/gridsize)*gridsize:gridsize:floor(stop/gridsize)*gridsize,
  plot([i i], [step/4 step*length(truepos)+step/2], 'k:') ;
  hold on
  h=text(i, step/8, sprintf('%i', i)) ;
  set(h, 'HorizontalAlignment', 'center') ;
  set(h, 'FontSize', 7) ;
end ;


if nargin<4 | isempty(colors)
  colors={} ;
  for i=1:500,
    colors{i,1}='k' ;
    colors{i,2}=[0 0 0] ;
  end ;
end ;
if nargin<3,
  labels={'Annotation', 'EST', 'Prediction'} ;
end ;


h=text(start, cnt+(length(truepos)+1)*step, name) ;
set(h, 'HorizontalAlignment', 'center') ;
set(h, 'FontSize', 10) ;
for k=1:length(truepos)
  
  cnt=cnt+step ;
  
  plot([start truepos{k}(1,1)], [cnt cnt], colors{k,1}) ;
  hold on
  plot([truepos{k}(end,2) stop], [cnt cnt], colors{k,1}) ;
  axis off
  h=text(start, cnt, [labels{k} ' ']) ;
%  set(h, 'VerticalAlignment', 'center') ;
  set(h, 'HorizontalAlignment', 'right') ;
  set(h, 'FontSize', 7) ;
  h=text(stop, cnt, ['  ' strands(k)]) ;
%  set(h, 'VerticalAlignment', 'center') ;
  set(h, 'HorizontalAlignment', 'left') ;
  set(h, 'FontSize', 7) ;
  axis([start stop 0 cnt+step]) ;
  
  for i=1:size(truepos{k},1)
    if size(truepos{k},2)==2 | truepos{k}(i,3)==1
      patch([truepos{k}(i,1) truepos{k}(i,1) truepos{k}(i,2) truepos{k}(i,2)],...
            [cnt-width, cnt+width, cnt+width, cnt-width], colors{k,1}) ;
    else
      patch([truepos{k}(i,1) truepos{k}(i,1) truepos{k}(i,2) truepos{k}(i,2)],...
            [cnt-width2, cnt+width2, cnt+width2, cnt-width2], colors{k,1}) ;
    end ;
    h=text(mean([truepos{k}(i,1) truepos{k}(i,2)]), cnt-3*width, sprintf('%i', truepos{k}(i,2)-truepos{k}(i,1))) ;
    set(h, 'HorizontalAlignment', 'center') ;
    set(h, 'FontSize', 7) ;
  end ;
  
  for i=1:size(truepos{k},1)-1
    h=line([truepos{k}(i,2) mean([truepos{k}(i,2) truepos{k}(i+1,1)]) truepos{k}(i+1,1)],...
         [cnt+width/2, cnt+2*width, cnt+width/2]) ;
    set(h,'color',colors{k,2}) ;
    h=text(mean([truepos{k}(i,2) truepos{k}(i+1,1)]), cnt+4*width, sprintf('%i', truepos{k}(i+1,1)-truepos{k}(i,2))) ;
    set(h, 'HorizontalAlignment', 'center') ;
    set(h, 'FontSize', 7) ;
  end ;
end ;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






cnt = cnt+3*step;

h=text(mean(start,stop),cnt,'Location of interest') ;
set(h, 'HorizontalAlignment', 'center') ;
set(h, 'FontSize', 10) ;

exon
strands
patch([exon(1) exon(1) exon(2) exon(2)],...
      [cnt-width, cnt+width, cnt+width, cnt-width], 'g');


cnt = cnt+step;

h=text(mean(start,stop),cnt,'Splice Graph') ;
set(h, 'HorizontalAlignment', 'center') ;
set(h, 'FontSize', 10) ;


exon_loc = zeros(1,stop-start);
exon_level = zeros(size(vertices,2));

for i=1:size(vertices,2)
  cur_vertex = vertices(:,i) - [start;start];
  exon_loc(cur_vertex(1):cur_vertex(2)) = ...
      exon_loc(cur_vertex(1):cur_vertex(2)) + 1;
  exon_loc(cur_vertex(1):cur_vertex(2)) = ...
      max(exon_loc(cur_vertex(1):cur_vertex(2)));
      
  level = max(exon_loc(cur_vertex(1):cur_vertex(2)));
  
  exon_level(i) = level;
  
  %disp([i,cur_vertex'+[start,start],level])

  
  
  patch([cur_vertex(1)+start cur_vertex(1)+start,...
	 cur_vertex(2)+start cur_vertex(2)+start],...
	[cnt+(level*step)-width, cnt+(level*step)+width,...
	 cnt+(level*step)+width, cnt+(level*step)-width], 'r');

  
  h=text(mean([vertices(1,i),vertices(2,i)]), cnt+(level*step)-3*width,...
	 sprintf('%i', i)) ;
  set(h, 'HorizontalAlignment', 'center') ;
  set(h, 'FontSize', 7) ;

  
end

axis([start stop 0 cnt+step*(max([exon_loc])+1)]) ;


intron_loc = zeros(1,stop-start);
if (size(edges,1) > 1)

  for i=1:size(vertices,2)
    for j=i+1:size(vertices,2)
      if edges(i,j)
	if (vertices(1,i) < vertices(1,j))
	  istart = vertices(2,i);
	  istop =  vertices(1,j);
	  level1 = exon_level(i);
	  level2 = exon_level(j);
	else
	  istart = vertices(2,j);
	  istop =  vertices(1,i);
	  level1 = exon_level(j);
	  level2 = exon_level(i);
	end
	
%	disp([i,j,istart,istop])
	
	
	
	cur_intron = [istart-start,istop-start];
	intron_loc(cur_intron(1):cur_intron(2)) = ...
	    intron_loc(cur_intron(1):cur_intron(2)) + 1;
	leveli = max(intron_loc(cur_intron(1):cur_intron(2)));
      
%      disp([i,j,istart,istop,exon_level(i),exon_level(j)])

      
%      line([istart, mean([istart,istop]),istop],...
%	   [cnt+(level1*step)+width/2,...
%	    mean([cnt+(level1*step)+2*width,cnt+(level2*step)+5*width]),...
%	    cnt+(level2*step)+width/2]);
      
        line([istart, mean([istart,istop]),istop],...
	     [cnt+(level1*step)+width/2,...
	      cnt+((leveli+1)*step*0.75),...
	      cnt+(level2*step)+width/2]);
      
      end
    end
  end
end


axis([start stop 0 cnt+step*(max([exon_loc,intron_loc])+1)]) ;
%axis([start stop cnt cnt+step*(max([exon_loc,intron_loc])+1)]) ;



return

