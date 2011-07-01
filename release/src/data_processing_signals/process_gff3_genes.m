function [genes, FAULTS_id,FAULTS_name] = process_gff3_genes(genes, genome_info, name_of_CDS,model) 
% genes = process_gff3_genes(genes, genome_info, name_of_CDS)  

if nargin<3
  name_of_CDS = 'CDS';
end

if nargin<4
  m.exon = 0;
  m.transexon = 1;
  m.utr5exon = 2;
  m.utr3exon = 4; 
  m.cds_exon = 100;
  if any(ismember(name_of_CDS, {'CDS','coding_exon'}))
    m.cds_exon = 3;
  end
  if any(ismember(name_of_CDS, {'exon'}))
    m.cds_exon = [m.cds_exon 0];
  end
else
 	if isfield(model.segments, 'exon')
		m.exon = model.segments.exon;
	else
		m.exon = -1e7;
	end
	if isfield(model.segments, 'trans_exon')
		m.transexon = model.segments.trans_exon;
	else
		m.transexon = -1e7; 
	end
	m.utr5exon = model.segments.utr5exon;
	m.utr3exon = model.segments.utr3exon;
	m.cds_exon = model.segments.cds_exon;
	if isfield(model.use, 'non_coding')&&model.use.non_coding
		m.nc_exon = model.segments.nc_exon;
	end
end

gff_contig_names = unique({genes.chr}) ;
unannotated_contigs = setdiff(upper(genome_info.contig_names), upper(gff_contig_names)) ;
if ~isempty(unannotated_contigs),
  fprintf(1, '\nThe genome information object contains %i/%i contigs not annotated in the GFF3 file.\n', length(unannotated_contigs), length(genome_info.contig_names)) ;
  if length(unannotated_contigs)<=10,
    for i=1:length(unannotated_contigs), 
      fprintf(1, '  %s\n', unannotated_contigs{i}) ;
    end;
  end ;
end ;
unknown_contigs = setdiff(upper(gff_contig_names), upper(genome_info.contig_names)) ;
if ~isempty(unknown_contigs),
  fprintf(2, '\nThe GFF3 file contains the following %i contigs not present in the genome information object:\n', length(unknown_contigs)) ;
  for i=1:length(unknown_contigs), 
    fprintf(2, '  %s\n', unknown_contigs{i}) ;
  end;
  error('Cannot continue.') ;
end ;


FAULTS_name.exons.overlapping = {};
FAULTS_name.exons.unknown_type = {};
FAULTS_name.exons.no_transcripts = [];
FAULTS_name.exons.no_exons = [];
FAULTS_name.exons.TAGT_cases = [];

FAULTS_name.cdsUTR5.ok = {};
FAULTS_name.cdsUTR5.TISdiff = {};
FAULTS_name.cdsUTR5.mult_overlap = {}; 
FAULTS_name.cdsUTR5.gap = {};
FAULTS_name.cdsUTR5.TSSdiff = {};
FAULTS_name.cdsUTR5.no_utr5 = {};

FAULTS_name.cdsUTR3.ok = {};
FAULTS_name.cdsUTR3.STOPdiff = {};
FAULTS_name.cdsUTR3.mult_overlap = {}; 
FAULTS_name.cdsUTR3.gap = {};
FAULTS_name.cdsUTR3.CLEAVEdiff = {};
FAULTS_name.cdsUTR3.no_utr3 = {};

FAULTS_name.exonUTR5.gap = {};
FAULTS_name.exonUTR5.TSSdiff = {};
FAULTS_name.exonUTR3.CLEAVEdiff = {};
FAULTS_name.exonUTR3.gap = {};
FAULTS_name.noCDS.noCDS = {};



FAULTS_id.exons.overlapping = [];
FAULTS_id.exons.unknown_type = [];
FAULTS_id.exons.no_transcripts = [];
FAULTS_id.exons.no_exons = [];
FAULTS_id.exons.TAGT_cases = [];

FAULTS_id.cdsUTR5.ok = [];
FAULTS_id.cdsUTR5.TISdiff = [];
FAULTS_id.cdsUTR5.mult_overlap = []; 
FAULTS_id.cdsUTR5.gap = [];
FAULTS_id.cdsUTR5.TSSdiff = [];
FAULTS_id.cdsUTR5.no_utr5 = [];

FAULTS_id.cdsUTR3.ok = [];
FAULTS_id.cdsUTR3.STOPdiff = [];
FAULTS_id.cdsUTR3.mult_overlap = []; 
FAULTS_id.cdsUTR3.gap = [];
FAULTS_id.cdsUTR3.CLEAVEdiff = [];
FAULTS_id.cdsUTR3.no_utr3 = [];

FAULTS_id.exonUTR5.gap = [];
FAULTS_id.exonUTR5.TSSdiff = [];
FAULTS_id.exonUTR3.CLEAVEdiff = [];
FAULTS_id.exonUTR3.gap = [];
FAULTS_id.noCDS.noCDS = [];

fprintf('processing genes\n')
for i=1:length(genes)
  if (mod(i,1000)==0) || i==length(genes), 
    fprintf('%i/%i\r',i,length(genes));
  end
  genes(i).id = i;
  if exist('genome_info','var') && ~isempty(genome_info) && (~isfield(genes,'chr_num') || isempty(genes(i).chr_num))
    genes(i).chr_num= strmatch(upper(genes(i).chr), upper(genome_info.contig_names), 'exact');
    if isempty(genes(i).chr_num)
      error('sequence ID in gff file (%s) does not match any sequence IDs in fasta file (%s)', upper(genes(i).chr), upper(genome_info.contig_names));
    end
  end
  genes(i).tis = {} ;
  genes(i).cdsStop = {} ;
  genes(i).tss = {} ;
  genes(i).cleave = {} ;
  if ~isfield(genes(i), 'cds_exons') || isempty(genes(i).cds_exons)
    genes(i).cds_exons = {} ;
  end ;
  if ~isfield(genes(i), 'utr5_exons') || isempty(genes(i).utr5_exons)
    genes(i).utr5_exons = {} ;
  end ;
  if ~isfield(genes(i), 'utr3_exons') || isempty(genes(i).utr3_exons)
    genes(i).utr3_exons = {} ;
  end ;
  gene = genes(i) ;
  if isempty(gene.exons) || (length(gene.exons)==1&&isempty(gene.exons{1}))
    FAULTS_id.exons.no_exons(end+1)=i;
    FAULTS_name.exons.no_exons{end+1}=gene.name;
    continue
  end
  if isempty(gene.transcripts) 
    FAULTS_id.exons.no_transcripts(end+1) = i;
    FAULTS_name.exons.no_transcripts{end+1} = gene.name;
    continue
  end
  
  for j=1:length(gene.transcripts)
    gene.exons_confirmed{j} =[]; 
    exons = gene.exons{j} ;
    if size(exons,2)<3, 
      %warning('no exon meta info, identifying as exons') ;
      exons(:,3)=m.exon ; 
    end ;
    [exons(:,1),idx] = sort(exons(:,1)) ;
    exons(:,2:end) = exons(idx,2:end) ;
    %assert(all((exons(:,2)-exons(:,1))>=0)) ;
    if ~all((exons(:,2)-exons(:,1))>=0),
      warning('exon_end<exon_start; gene invalid') ;
      continue ;
    end ;
      
    idx_trans = find(exons(:,3)==m.transexon); 
    if ~isempty(idx_trans)
      exons(idx_trans,:) = [];
      warning('transexon removed')
    end
    %remove redundant information 
    if all(ismember([m.exon,m.cds_exon],unique(exons(:,3))))
      idx_cds = find(exons(:,3)==m.cds_exon);
      idx_exon = find(exons(:,3)==m.exon);
      [tmp,tmp,idx2] = intersect(exons(idx_cds,1:2),exons(idx_exon,1:2),'rows');
      exons(idx_exon(idx2),:)=[];
    end
    %% overlapping exons are merged
    if ~all((exons(2:end,1)-exons(1:end-1,2))>=0)
      idx = find(exons(2:end,1)-exons(1:end-1,2)<0) ;
      if (exons(idx,3)==exons(idx+1,3))
        FAULTS_id.exons.overlapping(end+1) = i;
        FAULTS_name.exons.overlapping{end+1} = gene.name;
        % warning('overlapping exons found within one transcript\n')  
	if keyboard_allowed(),
          keyboard
        end ;
	exons(idx,2) =  exons(idx+1,2);
        exons(idx+1,:) =[];
        if ~all((exons(:,2)-exons(:,1))>=0)
	  error('corrupted exons found: length of exon is <= zero');
	end
	if ~all((exons(2:end,1)-exons(1:end-1,2))>=0)
	  error('overlapping exons found within a single transcript');
	end
      end
    end	
    
    all_exons = exons;
    cds_exons = exons(find(ismember(exons(:,3),m.cds_exon)),1:2) ;
    utr5_exons = exons(find(exons(:,3)==m.utr5exon), 1:2) ;
    utr3_exons = exons(find(exons(:,3)==m.utr3exon), 1:2) ;
    other = exons(find(~ismember(exons(:,3),[m.cds_exon m.utr5exon m.utr3exon m.exon])),:) ;
    exons = exons(find(exons(:,3)==m.exon),1:2) ;
    assert(size(all_exons,1)==size(cds_exons,1)+size(utr5_exons,1)+size(utr3_exons,1)+size(exons,1)+size(other,1))
    
    if size(other,1)>0
      FAULTS_name.exons.unknown_type{end+1} = gene.name;
      FAULTS_id.exons.unknown_type(end+1) = i;  
      %warning('unknown childs ignored')
    end
    tis = [];
    Stop = [];
    tss = [];
    cleave = [];
    if gene.strand=='+',
      
      if ~isempty(cds_exons)  
        tis = cds_exons(1,1) ;
        Stop = cds_exons(end,2) ;
        if ~isempty(utr5_exons)  %% 
          if utr5_exons(end,2)==cds_exons(1,1)
            %% all good
            FAULTS_id.cdsUTR5.ok(end+1) = i;
            FAULTS_name.cdsUTR5.ok{end+1} = gene.name;
          else
            FAULTS_name.cdsUTR5.TISdiff{end+1} = gene.name;
            FAULTS_id.cdsUTR5.TISdiff(end+1) = i; 
            % warning(sprintf('last UTR-5 exon end (%i) is not equal to the start of the first coding exon (%i) (gene %s, contig %s)\nignoring 5''UTR\n', utr5_exons(end,2), cds_exons(1,1), gene.name, gene.chr));
            utr5_exons = [] ;
          end
	elseif ~isempty(exons) && any(exons(:,1)<tis) 
          idx = find(exons(:,1)<tis);
          utr5_exons = exons(idx,:);
          if sum(utr5_exons(:,2)>tis)>1
            FAULTS_name.cdsUTR5.mult_overlap{end+1} = gene.name;
            FAULTS_id.cdsUTR5.mult_overlap(end+1) = i; 
            % warning(sprintf('more than one exons overlaps start of the first coding exon (gene %s, contig %s) \nignoring 5''UTR \n', gene.name, gene.chr));
            utr5_exons = [];
          elseif ~any(utr5_exons(:,2)>tis)
            utr5_exons(end,2) = cds_exons(1,1);
            FAULTS_name.cdsUTR5.gap{end+1} = gene.name;
            FAULTS_id.cdsUTR5.gap(end+1) = i; 
            % warning(sprintf('extending 5'' UTR to start of CDS exon\n', gene.name, gene.chr));
          elseif utr5_exons(end,2)~=cds_exons(1,2)
            FAULTS_name.cdsUTR5.TSSdiff{end+1} = gene.name;
            FAULTS_id.cdsUTR5.TSSdiff(end+1) = i; 
            % warning(sprintf('TIS overlapping exon has different end than first cds_exon\n', gene.name, gene.chr));
            utr5_exons(end,2) = cds_exons(1,1);
          else
            FAULTS_id.cdsUTR5.ok(end+1) = i;
            FAULTS_name.cdsUTR5.ok{end+1} = gene.name;
            utr5_exons(end,2) = cds_exons(1,1);
          end
        else
          FAULTS_name.cdsUTR5.no_utr5{end+1} = gene.name;
          FAULTS_id.cdsUTR5.no_utr5(end+1) = i; 
          % warning(sprintf('cannot find 5'' UTR \n', gene.name, gene.chr));
        end ;
        
        if ~isempty(utr3_exons)  %% 
          if utr3_exons(1,1)==cds_exons(end,2)
            %% all good
            FAULTS_id.cdsUTR3.ok(end+1) = i;
            FAULTS_name.cdsUTR3.ok{end+1} = gene.name;
          else
            FAULTS_name.cdsUTR3.STOPdiff{end+1} = gene.name;
            FAULTS_id.cdsUTR3.STOPdiff(end+1) = i; 
            % warning(sprintf('first UTR-3 exon start (%i) is not equal to the end of the first coding exon (%i) (gene %s, contig %s)\nignoring 3''UTR\n', utr3_exons(end,2), cds_exons(1,1), gene.name, gene.chr));
            utr3_exons = [] ;
          end
        elseif ~isempty(exons) && any(exons(:,2)>Stop) 
          idx = find(exons(:,2)>Stop);
          utr3_exons = exons(idx,:);
          if sum(utr3_exons(:,1)<Stop)>1
            FAULTS_name.cdsUTR3.mult_overlap{end+1} = gene.name;
            FAULTS_id.cdsUTR3.mult_overlap(end+1) = i; 
            % warning(sprintf('more than one exons overlaps end of the first coding exon (gene %s, contig %s) \nignoring 3''UTR \n', gene.name, gene.chr));
            utr3_exons = [];
          elseif ~any(utr3_exons(:,1)<Stop)
            utr3_exons(1,1) = cds_exons(end,2);
            FAULTS_name.cdsUTR3.gap{end+1} = gene.name;
            FAULTS_id.cdsUTR3.gap(end+1) = i; 
            % warning(sprintf('extending 3'' UTR to end of CDS exon\n', gene.name, gene.chr));
          elseif utr3_exons(1,1)~=cds_exons(end,1)
            FAULTS_name.cdsUTR3.CLEAVEdiff{end+1} = gene.name;
            FAULTS_id.cdsUTR3.CLEAVEdiff(end+1) = i; 
            % warning(sprintf('TIS overlapping exon has different end than first cds_exon\n', gene.name, gene.chr));
            utr3_exons(1,1) = cds_exons(end,2);
          else
            FAULTS_id.cdsUTR3.ok(end+1) = i;
            FAULTS_name.cdsUTR3.ok{end+1} = gene.name;
            utr3_exons(1,1) = cds_exons(end,2);
          end
        else
          FAULTS_name.cdsUTR3.no_utr3{end+1} = gene.name;
          FAULTS_id.cdsUTR3.no_utr3(end+1) = i; 
          % warning(sprintf('cannot find 3'' UTR \n', gene.name, gene.chr));
        end ; % isempty(utr3_exons)
        
      % start no_cds_exons given
      elseif ~isempty(utr5_exons) & ~isempty(exons) & ~isempty(utr3_exons)
        % inferring cds_exons
        tis = utr5_exons(end,2) ;
        Stop = utr3_exons(1,1) ;
        idx = find(exons(:,2)>tis&exons(:,1)<Stop);
        cds_exons = exons(idx,:);
        if cds_exons(1,1)<tis 
          if cds_exons(1,1)~= utr5_exons(end,1) ;
            FAULTS_name.exonUTR5.TSSdiff{end+1} = gene.name;
            FAULTS_id.exonUTR5.TSSdiff(end+1) = i; 
            % warning(['first potential coding exon has different start than ' ...
            %          'last 5''UTR exon']) 
          end
          cds_exons(1,1)=tis;
        else
          FAULTS_name.exonUTR5.gap{end+1} = gene.name;
          FAULTS_id.exonUTR5.gap(end+1) = i; 
          % warning(['first potential coding exon starts downstream of tis']) 
          cds_exons(1,1)=tis;
        end
        if cds_exons(end,2)>Stop 
          if cds_exons(end,2)~= utr3_exons(1,2) ;
            FAULTS_name.exonUTR3.CLEAVEdiff{end+1} = gene.name;
            FAULTS_id.exonUTR3.CLEAVEdiff(end+1) = i; 
            % warning(['last potential coding exon has different end than ' ...
            %          'first 3''UTR exon']) 
          end
          cds_exons(end,2)=Stop;
        else
          FAULTS_name.exonUTR3.gap{end+1} = gene.name;
          FAULTS_id.exonUTR3.gap(end+1) = i; 
          % warning(['last potential coding exon ends upstream of cdsStop']) 
          cds_exons(end,2)=Stop;
        end 
      else
        tis = [] ;
        Stop = [] ;
        FAULTS_name.noCDS.noCDS{end+1} = gene.name;
        FAULTS_id.noCDS.noCDS(end+1) = i; 
        % warning(sprintf('cannot infer UTRs and coding regions \n', gene.name, gene.chr));
      end
      
      if ~isempty(utr3_exons) || ~isempty(utr5_exons) || ~isempty(cds_exons),
        Stop=Stop-3 ;
        
        if size(utr3_exons,1)>1 && utr3_exons(1,2)-utr3_exons(1,1)<3
          utr3_exons = [utr3_exons(1,1) utr3_exons(end,2)];
          FAULTS_name.exons.TAGT_cases{end+1} = gene.name;
          FAULTS_id.exons.TAGT_cases(end+1) = i;
          % warning('TAGT case: invalid donor!')
        end
        % determine phase of cds exon
        cds_exons(1,3) = 0;
        cds_length = 0;
        for e=2:size(cds_exons,1)
          cds_length = cds_length + cds_exons(e-1,2)-cds_exons(e-1,1);
          phase = mod(cds_length,3);
          switch phase
           case 0
            cds_exons(e,3) = phase;
           case 1
            cds_exons(e,3) = 2;
           case 2
            cds_exons(e,3) = 1;
          end
        end
        exons = utr5_exons;
        if ~isempty(utr5_exons)
          tss =  utr5_exons(1,1);
        end ;
        if ~isempty(cds_exons)
          if ~isempty(exons)
            exons(end,2) = cds_exons(1,2) ;
            exons = [exons; cds_exons(2:end,1:2)] ;
          else
            exons = cds_exons(:,1:2) ;
          end ;
        end
        if ~isempty(utr3_exons)
          cleave =  utr3_exons(end,2);
          if ~isempty(exons)
            exons(end,2) = utr3_exons(1,2) ;
            exons = [exons; utr3_exons(2:end,:)] ;
          else
            exons = utr3_exons(:,1:2) ;
          end ;
        end
      else
        if ~isempty(exons)
          tss = exons(1,1) ;
          cleave = exons(end,2) ;
        end ;
      end 
      
    else % gene.strand=='-',
      if ~isempty(cds_exons)  
        tis = cds_exons(end,2) ;
        Stop = cds_exons(1,1) ;
        if ~isempty(utr5_exons)  %% 
          if utr5_exons(1,1)==cds_exons(end,2)
            %% all good
            FAULTS_id.cdsUTR5.ok(end+1) = i;
            FAULTS_name.cdsUTR5.ok{end+1} = gene.name;
          else
            FAULTS_name.cdsUTR5.TISdiff{end+1} = gene.name;
            FAULTS_id.cdsUTR5.TISdiff(end+1) = i; 
            % warning(sprintf('last UTR-5 exon end (%i) is not equal to the start of the first coding exon (%i) (gene %s, contig %s)\nignoring 5''UTR\n', utr5_exons(end,2), cds_exons(1,1), gene.name, gene.chr));
            utr5_exons = [] ;
          end
	elseif ~isempty(exons) && any(exons(:,2)>tis) 
          idx = find(exons(:,2)>tis);
          utr5_exons = exons(idx,:);
          if sum(utr5_exons(:,1)<tis)>1
            FAULTS_name.cdsUTR5.mult_overlap{end+1} = gene.name;
            FAULTS_id.cdsUTR5.mult_overlap(end+1) = i; 
            % warning(sprintf('more than one exons overlaps start of the first coding exon (gene %s, contig %s) \nignoring 5''UTR \n', gene.name, gene.chr));
            utr5_exons = [];
          elseif ~any(utr5_exons(:,1)<tis)
            FAULTS_name.cdsUTR5.gap{end+1} = gene.name;
            FAULTS_id.cdsUTR5.gap(end+1) = i; 
            utr5_exons(1,1) = cds_exons(end,2);
            % warning(sprintf('extending 5'' UTR to start of CDS exon\n', gene.name, gene.chr));
          elseif utr5_exons(1,1)~=cds_exons(end,1)
            FAULTS_name.cdsUTR5.TSSdiff{end+1} = gene.name;
            FAULTS_id.cdsUTR5.TSSdiff(end+1) = i; 
            % warning(sprintf('TIS overlapping exon has different end than first cds_exon\n', gene.name, gene.chr));
            utr5_exons(1,1) = cds_exons(end,2);
          else
            utr5_exons(1,1) = cds_exons(end,2);
          end
        else
          FAULTS_name.cdsUTR5.no_utr5{end+1} = gene.name;
          FAULTS_id.cdsUTR5.no_utr5(end+1) = i; 
          % warning(sprintf('cannot find 5'' UTR \n', gene.name, gene.chr));
        end ;
        
        if ~isempty(utr3_exons)  %% 
          if utr3_exons(end,2)==cds_exons(1,1)
            %% all good
            FAULTS_id.cdsUTR3.ok(end+1) = i;
            FAULTS_name.cdsUTR3.ok{end+1} = gene.name;
          else
            FAULTS_name.cdsUTR3.STOPdiff{end+1} = gene.name;
            FAULTS_id.cdsUTR3.STOPdiff(end+1) = i; 
            %  warning(sprintf('first UTR-3 exon start (%i) is not equal to the end of the first coding exon (%i) (gene %s, contig %s)\nignoring 5''UTR\n', utr3_exons(end,2), cds_exons(1,1), gene.name, gene.chr));
            utr3_exons = [] ;
          end
        elseif ~isempty(exons) && any(exons(:,1)<Stop) 
          idx = find(exons(:,1)<Stop);
          utr3_exons = exons(idx,:);
          if sum(utr3_exons(:,2)>Stop)>1
            FAULTS_name.cdsUTR3.mult_overlap{end+1} = gene.name;
            FAULTS_id.cdsUTR3.mult_overlap(end+1) = i; 
            % warning(sprintf('more than one exons overlaps end of the first coding exon (gene %s, contig %s) \nignoring 5''UTR \n', gene.name, gene.chr));
            utr3_exons = [];
          elseif ~any(utr3_exons(:,2)>Stop)
            FAULTS_name.cdsUTR3.gap{end+1} = gene.name;
            FAULTS_id.cdsUTR3.gap(end+1) = i; 
            utr3_exons(end,2) = cds_exons(1,1);
            % warning(sprintf('extending 3'' UTR to end of CDS exon\n', gene.name, gene.chr));
          elseif utr3_exons(end,2)~=cds_exons(1,2)
            FAULTS_name.cdsUTR3.CLEAVEdiff{end+1} = gene.name;
            FAULTS_id.cdsUTR3.CLEAVEdiff(end+1) = i; 
            % warning(sprintf('TIS overlapping exon has different end than first cds_exon\n', gene.name, gene.chr));
            utr3_exons(end,2) = cds_exons(1,1);
          else
            utr3_exons(end,2) = cds_exons(1,1);
          end
        else
          FAULTS_name.cdsUTR3.no_utr3{end+1} = gene.name;
          FAULTS_id.cdsUTR3.no_utr3(end+1) = i; 
          % warning(sprintf('cannot find 3'' UTR \n', gene.name, gene.chr));
        end ;

      % start no_cds_exons given
      elseif ~isempty(utr5_exons) & ~isempty(exons) & ~isempty(utr3_exons)
        % inferring cds_exons
        tis = utr5_exons(1,1) ;
        Stop = utr3_exons(end,2) ;
        idx = find(exons(:,1)<tis&exons(:,2)>Stop);
        cds_exons = exons(idx,:);
        if cds_exons(end,2)>tis 
          if cds_exons(end,2)~= utr5_exons(1,2) ;
            FAULTS_name.exonUTR5.TSSdiff{end+1} = gene.name;
            FAULTS_id.exonUTR5.TSSdiff(end+1) = i; 
            % warning(['first potential coding exon has different start than ' ...
            %          'last 5''UTR exon']) 
          end
          cds_exons(end,2)=tis;
        else
          FAULTS_name.exonUTR5.gap{end+1} = gene.name;
          FAULTS_id.exonUTR5.gap(end+1) = i; 
          % warning(['first potential coding exon starts downstream of tis']) 
          cds_exons(end,2)=tis;
        end
        if cds_exons(1,1)<Stop 
          if cds_exons(1,1)~= utr3_exons(end,1) ;
            FAULTS_name.exonUTR3.CLEAVEdiff{end+1} = gene.name;
            FAULTS_id.exonUTR3.CLEAVEdiff(end+1) = i; 
            % warning(['last potential coding exon has different end than ' ...
            %          'first 3''UTR exon']) 
          end
          cds_exons(1,1)=Stop;
        else
          FAULTS_name.exonUTR3.gap{end+1} = gene.name;
          FAULTS_id.exonUTR3.gap(end+1) = i; 
          % warning(['last potential coding exon ends upstream of Stop']) 
          cds_exons(1,1)=Stop;
        end 
      else
        tis = [] ;
        Stop = [] ;
        FAULTS_name.noCDS.noCDS{end+1} = gene.name;
        FAULTS_id.noCDS.noCDS(end+1) = i; 
        % warning(sprintf('cannot infer UTRs and coding regions \n', gene.name, gene.chr));
      end
      
      if ~isempty(utr3_exons) || ~isempty(utr5_exons) || ~isempty(cds_exons),
        Stop=Stop+3 ;
        
        if size(utr3_exons,1)>1 && utr3_exons(end,2)-utr3_exons(end,1)<3
          utr3_exons = [utr3_exons(1,1) utr3_exons(end,2)];
          FAULTS_name.exons.TAGT_cases{end+1} = gene.name;
          FAULTS_id.exons.TAGT_cases(end+1) = i;
          % warning('TAGT case: invalid donor!')
        end
        
        % determine phase of cds exon
        if ~isempty(cds_exons)
          cds_exons(end,3) = 0;
          cds_length = 0;
          for e=size(cds_exons,1)-1:-1:1
            cds_length = cds_length + cds_exons(e+1,2)-cds_exons(e+1,1);
            phase = mod(cds_length,3);
            switch phase
             case 0
              cds_exons(e,3) = phase;
             case 1
              cds_exons(e,3) = 2;
             case 2
              cds_exons(e,3) = 1;
            end
          end
        end    
        exons = utr3_exons;
        if ~isempty(utr3_exons)
          cleave = utr3_exons(1,1);
        end 
        if ~isempty(cds_exons)
          if ~isempty(exons),
            exons(end,2) = cds_exons(1,2) ;
            exons = [exons; cds_exons(2:end,1:2)] ;
          else
            exons = cds_exons(:,1:2) ;
          end ;
        end
        if ~isempty(utr5_exons)
          tss = utr5_exons(end,2);
          if ~isempty(exons)
            exons(end,2) = utr5_exons(1,2) ;
            exons = [exons; utr5_exons(2:end,:)] ;
          else
            exons = utr5_exons(:,1:2) ;
          end 
        end
      else
        if ~isempty(exons),
          tss = exons(end,2) ;
          cleave = exons(1,1) ;
        end ;
      end 
      
    end ; %strand==-'

    gene.exons{j} = exons ;
    gene.cds_exons{j} = cds_exons;
    gene.utr5_exons{j} = utr5_exons;
    gene.utr3_exons{j} = utr3_exons;
    
    gene.tis{j}=tis ;
    gene.cdsStop{j}=Stop ;
    gene.tss{j} = tss;
    gene.cleave{j} = cleave;
    
  end ; %%loop over transcripts
  
  gene.start = min(gene.exons{1}(:)) ;
  gene.stop = max(gene.exons{1}(:)) ;
  for j=2:length(gene.exons)
    gene.start = min(gene.start, min(gene.exons{j}(:))) ;
    gene.stop = max(gene.stop, max(gene.exons{j}(:))) ;
  end ;
  
  genes(i)=gene ;
end ;%%loop over genes



F1=fieldnames(FAULTS_name);
for f=1:length(F1)
  F=fieldnames(FAULTS_name.(F1{f}));
  for i=1:length(F)
    [tmp,idx] = unique(FAULTS_id.(F1{f}).(F{i}));
    FAULTS_name.(F1{f}).(F{i}) = FAULTS_name.(F1{f}).(F{i})(idx);
    FAULTS_id.(F1{f}).(F{i})= FAULTS_id.(F1{f}).(F{i})(idx);
  end
end

fprintf('\n\n');
fprintf('number of genes with overlapping exons:\t\t\t\t%i \n',length(unique(FAULTS_id.exons.overlapping)));
fprintf('number of genes with exons of unknown type:\t\t\t%i \n',length(unique(FAULTS_id.exons.unknown_type)));
fprintf('number of genes with invalid donor (TAGT):\t\t\t%i \n' ,length(unique(FAULTS_name.exons.TAGT_cases)));
        
fprintf('number of genes with cds and UTR5 derived TIS different:\t%i \n',length(unique(FAULTS_id.cdsUTR5.TISdiff)));
fprintf('number of genes with multiple UTR5s overlapping TIS: \t\t%i \n',length(unique(FAULTS_id.cdsUTR5.mult_overlap)));
fprintf('number of genes with gap between 5'' UTR and CDS: \t\t%i \n',length(unique(FAULTS_id.cdsUTR5.gap)));
fprintf('number of genes with different end of first CDS exon:\t\t%i \n',length(unique(FAULTS_id.cdsUTR5.TSSdiff)));
fprintf('number of genes with coding region but no 5'' UTR:\t\t%i \n',length(unique(FAULTS_id.cdsUTR5.no_utr5)));

fprintf('number of genes with cds and UTR3 derived STOP different:\t%i \n',length(unique(FAULTS_id.cdsUTR3.STOPdiff)));
fprintf('number of genes with multiple UTR3s overlapping STOP:\t\t%i \n',length(unique(FAULTS_id.cdsUTR3.mult_overlap)));
fprintf('number of genes with gap between 3'' UTR and CDS:\t\t%i \n',length(unique(FAULTS_id.cdsUTR3.gap)));
fprintf('number of genes with different start of last CDS exon:\t\t%i \n',length(unique(FAULTS_id.cdsUTR3.CLEAVEdiff)));
fprintf('number of genes with coding region but no 3'' UTR:\t\t%i \n',length(unique(FAULTS_id.cdsUTR3.no_utr3)));

fprintf('number of genes with exon and UTR5 derived TSS different:\t%i \n',length(unique(FAULTS_id.exonUTR5.TSSdiff)));
fprintf('number of genes with gap between 5'' UTR and exon:\t\t%i \n',length(unique(FAULTS_id.exonUTR5.gap)));
fprintf('number of genes with exon and UTR3 derived CLEAVE different:\t%i \n',length(unique(FAULTS_id.exonUTR3.CLEAVEdiff)));
fprintf('number of genes with gap between 3'' UTR and exon:\t\t%i \n',length(unique(FAULTS_id.exonUTR3.gap)));
fprintf('number of genes without coding or UTR information:\t\t%i \n',length(unique(FAULTS_id.noCDS.noCDS)));


