function [max_exon_level,max_intron_level,exon_level,intron_level,level] = detectsplicegraph(gene)

  [max_exon_level,max_intron_level,exon_level,intron_level,level] = countlevel(gene.exons,gene.name,...
                                                    gene.transcripts,[], ...
                                                    gene.splicegraph{1},...
                                                    gene.splicegraph{2});
  return


function [max_level,max_leveli,exon_level,intron_level,level] = countlevel(truepos,name,labels,colors,...
                                                    vertices,edges)
  
  
  start=inf ;
  stop=-inf ;
  for k=1:length(truepos)
    start=min(start,truepos{k}(1,1)-50) ;
    stop=max(stop,truepos{k}(end,2)+50) ;
  end ;
  
  
  max_level = 0;
  max_leveli = 0;

  splice = unique(vertices(:)) ;
  
  exon_loc = zeros(1,stop-start);
  exon_level = zeros(3,size(vertices,2));
  for i=1:size(vertices,2)
    cur_vertex = vertices(:,i) - [start;start];
    exon_loc(cur_vertex(1):cur_vertex(2)-1) = ...
        exon_loc(cur_vertex(1):cur_vertex(2)-1) + 1;
    %exon_loc(cur_vertex(1):cur_vertex(2)-1) = ...
    %    max(exon_loc(cur_vertex(1):cur_vertex(2)-1));
    
    level = max(exon_loc(cur_vertex(1):cur_vertex(2)-1));
    if isempty(level), level=0 ; end ;
    max_level = max(level,max_level);
    exon_level(1,i) = vertices(1,i) ;  
    exon_level(2,i) = vertices(2,i) ;  
    exon_level(3,i) = level;  
  end
  
  intron_loc = zeros(1,stop-start);
  intron_level = zeros(3,0) ;
  if (size(edges,1) > 1)
    for i=1:size(vertices,2)
      for j=i+1:size(vertices,2)
        if edges(i,j)
          if (vertices(1,i) < vertices(1,j))
            istart = vertices(2,i);
            istop =  vertices(1,j);
            %level1 = exon_level(i);
            %level2 = exon_level(j);
          else
            istart = vertices(2,j);
            istop =  vertices(1,i);
            %level1 = exon_level(j);
            %level2 = exon_level(i);
          end
          
          cur_intron = [istart-start,istop-start];
          intron_loc(cur_intron(1):cur_intron(2)-1) = ...
              intron_loc(cur_intron(1):cur_intron(2)-1) + 1;
          intron_level(3,end+1) = max(intron_loc(cur_intron(1):cur_intron(2)-1));
          intron_level(1,end) = istart ;
          intron_level(2,end) = istop ;

          max_leveli = max(intron_level(3,end),max_leveli);
        end
      end
    end
  end
  
  splice_loc = exon_loc+intron_loc ;
  level=zeros(5,length(splice)-1) ;
  for i=1:length(splice)-1,
    level(:,i)=[splice(i);
                splice(i+1);
                max(splice_loc([splice(i)-start:splice(i+1)-start-1]));
                max(exon_loc([splice(i)-start:splice(i+1)-start-1]));
                max(intron_loc([splice(i)-start:splice(i+1)-start-1]));
               ] ;
  end ;

return

