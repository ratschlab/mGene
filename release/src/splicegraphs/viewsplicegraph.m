function viewsplicegraph(gene,flag,param)

% viewsplicegraph(gene,flag)
% set flag = 1 to view only the graph
% set flag = 2 to view only the transcripts
% default behaviour views both.
%
% param is an optional structure with fields and default values:
%   text_on      1
%   FontSize     8
%   FontWeight   'normal'
%   LineWidth    1
%   grid_on      1
%   colors       {[0 0 0]   transcripts exon 
%                 [0 0 0]   transcripts intron
%                 [1 0 0]   summary exon
%                 [0 0 0]}  summary intron
%
% modified from cheng's version by GS 17/08/06
% => intron retention represented with both exons on the same level    
% => param added 
  
if nargin<3,param=struct; end ;
if nargin<2,flag=0 ; end ;

if ~isfield(gene,'strands')
  gene.strands = gene.strand;
end

if flag==0,
  viewgraph(gene.exons,gene.name, gene.transcripts, gene.strands, ...
            gene.splicegraph{1},gene.splicegraph{2},gene.chr,param);
elseif flag==1,
  viewgraph({},gene.name, {}, gene.strands, ...
            gene.splicegraph{1},gene.splicegraph{2},gene.chr,param);
elseif flag==2,
  viewgraph(gene.exons,gene.name,gene.transcripts,gene.strands, ...
            [],[],gene.chr,param); 
end ;


return




function viewgraph(truepos, name, labels,strands,vertices,edges,chr,param)

%%%%% SET PARAMETERS
if length(strands)<length(labels)
  strands = repmat(strands(1),1,length(labels)) ;
end
  
width=0.075 ; width2=0.03; step=0.4 ;

if ~exist('param')
  param = struct;
end

if ~isfield(param,'FontSize')
  param.FontSize = 8;
end
if ~isfield(param,'text_on')
  param.text_on = 1;
end
if ~isfield(param,'FontWeight')
  param.FontWeight = 'normal';
end
if ~isfield(param,'LineWidth')
  param.LineWidth = 1;
end 
if ~isfield(param,'grid_on')
  param.grid_on = 1;
end
if ~isfield(param,'colors')
  param.colors={} ;
  for i=1:500,
    param.colors{i,1}=[0 0.3 0] ;% transcripts exon 
    param.colors{i,2}=[0 0.3 0] ;% transcripts intron
    param.colors{i,3}=[1 0 0] ;% summary exon
    param.colors{i,4}=[0 0 0] ;% summary intron
  end 
end


start=inf ;
stop=-inf ;
for k=1:length(truepos)
  start=min(start,truepos{k}(1,1)-50) ;
  stop=max(stop,truepos{k}(end,2)+50) ;
end ;

for k=1:size(vertices,1)
  start=min(start,min(vertices(1,:)-50)) ;
  stop=max(stop,max(vertices(2,:)+50)) ;
end ;

hold on
axis off
cnt = 0;
if length(truepos)>0,   
  if nargin<3,
    labels={'Annotation', 'EST', 'Prediction'} ;
  end ;

  for k=1:length(truepos)
    cnt=cnt+step ;
    
%%%%%% DRAW LINE FROM GENE START TO START OF FIRST EXON 
    h = plot([start truepos{k}(1,1)], [cnt cnt], 'k') ;
    set(h, 'LineWidth',0.5) ;
    hold on
%%%%%% DRAW LINE FROM END OF LAST EXON TO GENE END 
    h = plot([truepos{k}(end,2) stop], [cnt cnt],'k') ;
    set(h, 'LineWidth',0.5) ;
    axis off
%%%%%% ADD TRANSCRIPT NAME 
    if param.text_on,
      h = text(start, cnt, [labels{k} ' ']) ;
      set(h, 'HorizontalAlignment', 'right') ;
      set(h, 'FontSize',param.FontSize ) ;
      set(h, 'FontWeight',param.FontWeight) ;
    end ;
%%%%%% ADD STRAND 
    if param.text_on,
      h = text(stop, cnt, ['  ' strands(k)]) ;
      set(h, 'HorizontalAlignment', 'left') ;
      set(h, 'FontSize',param.FontSize ) ;
      set(h, 'FontWeight',param.FontWeight) ;
    end ;
    axis([start stop 0 cnt+step]) ;
    
%%%%%% DRAW EXONS 
    for i=1:size(truepos{k},1)
      if size(truepos{k},2)==2 | truepos{k}(i,3)==1
        patch([truepos{k}(i,1) truepos{k}(i,1) truepos{k}(i,2) truepos{k}(i,2)],...
              [cnt-width, cnt+width, cnt+width, cnt-width],param.colors{mod(k,499)+1,1}) ;
      else
        patch([truepos{k}(i,1) truepos{k}(i,1) truepos{k}(i,2) truepos{k}(i,2)],...
              [cnt-width2, cnt+width2, cnt+width2, cnt-width2],param.colors{mod(k,499)+1,1}) ;
      end ;
%%%%%% ADD EXON LENGTH
      if i<size(truepos{k},1) && truepos{k}(i,2)==truepos{k}(i+1,1),
        if param.text_on,
          h = text(mean([truepos{k}(i,1) truepos{k}(i+1,2)]), cnt-2*width, sprintf('%i', truepos{k}(i+1,2)-truepos{k}(i,1))) ;
          set(h, 'HorizontalAlignment', 'center') ;
          set(h, 'FontSize',param.FontSize) ;
          set(h, 'FontWeight',param.FontWeight) ;
        end ;
      elseif i>1 && truepos{k}(i-1,2)==truepos{k}(i,1),
      else
        if param.text_on,
          h = text(mean([truepos{k}(i,1) truepos{k}(i,2)]), cnt-2*width, sprintf('%i', truepos{k}(i,2)-truepos{k}(i,1))) ;
          set(h, 'HorizontalAlignment', 'center') ;
          set(h, 'FontSize',param.FontSize) ;
          set(h, 'FontWeight',param.FontWeight) ;
        end ;
      end ;
    end ;
%%%%%% DRAW INTRONS 
    for i=1:size(truepos{k},1)-1
      if truepos{k}(i,2)==truepos{k}(i+1,1), continue ; end ;
      h = line([truepos{k}(i,2) mean([truepos{k}(i,2) truepos{k}(i+1,1)]) truepos{k}(i+1,1)],...
             [cnt+width/2, cnt+2*width, cnt+width/2]) ;
      set(h,'color',param.colors{mod(k,499)+1,2}) ;
      set(h, 'LineWidth',param.LineWidth) ;
%%%%%% ADD INTRON LENGTH
      if param.text_on,
        h = text(mean([truepos{k}(i,2) truepos{k}(i+1,1)]), cnt+0*width, sprintf('%i', truepos{k}(i+1,1)-truepos{k}(i,2))) ;
        set(h, 'HorizontalAlignment', 'center') ;
        set(h, 'FontSize',param.FontSize) ;
        set(h, 'FontWeight',param.FontWeight) ;
      end ;
    end ;
  end ;
end ;
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%%%%%% ADD SUMMARY SPLICE GRAPH  


%%%%%% DRAW EXONS 

if ~isempty(truepos) && ~isempty(vertices)
  cnt = cnt+1*step;
end
exon_num = zeros(1,stop-start);
exon_loc = zeros(1,stop-start);
exon_level = zeros(size(vertices,2));
% keyboard
for i=1:size(vertices,2)
  cur_vertex = vertices(:,i) - [start;start];
  exon_num(cur_vertex(1):cur_vertex(2)) = ...
      exon_num(cur_vertex(1):cur_vertex(2)) + 1;
  if all(exon_num<2)
    exon_loc = exon_num;
    level = 1;
  elseif max(exon_num)>size(exon_loc,1)
    exon_loc(end+1,:) = 0; 
    exon_loc(end,cur_vertex(1):cur_vertex(2)) = 1; 
    level = size(exon_loc,1);
  elseif max(exon_num)<=size(exon_loc,1)
    idx = min(find(all(exon_loc(:,cur_vertex(1):cur_vertex(2))==0,2))); 
    exon_loc(idx,cur_vertex(1):cur_vertex(2)) = 1;
    level = idx;
  end
 
  exon_level(i) = level;
  
  patch([cur_vertex(1)+start cur_vertex(1)+start,...
	 cur_vertex(2)+start cur_vertex(2)+start],...
	[cnt+(level*step)-width, cnt+(level*step)+width,...
	 cnt+(level*step)+width, cnt+(level*step)-width],param.colors{mod(k,499)+1,3} );

  if param.text_on,
    
    h=text(mean([vertices(1,i),vertices(2,i)]), cnt+(level*step)-2*width,...
           sprintf('%i', i)) ;
    set(h, 'HorizontalAlignment', 'center') ;
    set(h, 'FontSize',param.FontSize) ;
    set(h, 'FontWeight',param.FontWeight) ;
  end ;
end

%%%%%% DRAW INTRONS 
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

	cur_intron = [istart-start,istop-start];
	intron_loc(cur_intron(1):cur_intron(2)) = ...
	    intron_loc(cur_intron(1):cur_intron(2)) + 1;
        leveli = (level1+level2).*.6;
        h = line([istart, mean([istart,istop]),istop],...
	     [cnt+(level1*step)+width/2,...
	      cnt+((leveli+1)*step*0.75),...
	      cnt+(level2*step)+width/2]);
        set(h, 'LineWidth',param.LineWidth) ;
        set(h, 'Color',param.colors{mod(k,499)+1,4}) ;
      end
    end
  end
end
if ~isempty(exon_level)
  cnt = cnt+max(max(exon_level))*step;
end
%%%%% DRAW GRID
if param.grid_on 
  gridsize = round((stop-start)/5);
  for i=ceil(start/gridsize)*gridsize:gridsize:floor(stop/gridsize)* ...
        gridsize,
    plot([i i], [step/4 cnt+1*step], 'k:') ;
    % plot([i i], [step/4 step*length(truepos)+step/2], 'k:') ;
    hold on
    if param.text_on,      
      h=text(i, step/8, sprintf('%i', i)) ;
      set(h, 'HorizontalAlignment', 'center') ;
      set(h, 'FontSize',param.FontSize) ;
      set(h, 'FontWeight',param.FontWeight) ;
    end ;
  end ;
end
  
%%%%%% ADD TITLE 
%h = text(mean(start,stop),cnt+1.5*step, sprintf('Scaffold %s %c',chr(10:end),strands(1))) ;
if param.text_on,
  h = text(mean(start,stop),cnt+1.5*step, sprintf('%s\nChromosome %s %c',name,chr,strands(1))) ;
  set(h, 'HorizontalAlignment', 'left') ;
  set(h, 'FontSize',param.FontSize ) ;
  set(h, 'FontWeight',param.FontWeight) ;
end ;

% keyboard

axis([start stop 0 cnt+2*step]) ;
hold off
% axis([start stop 0 cnt+step*max((max([exon_loc,intron_loc])+0.5),3)]) ;


return

