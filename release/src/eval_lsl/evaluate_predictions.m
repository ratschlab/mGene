function [res,R,map] = evaluate_predictions(regions,genes_anno,genes_pred,genome_info,run_local,fn_anno,fn_pred,fn_output);

% [res,map] = evaluate_predictions(regions,genes_anno,genes_pred,genome_info,view_on);

% Computes sensitivity and specivficity on (cds-) nucleotide, (cds-)
% exon, (cds-) transcript and gene level for a set of predicted genes,
% relative to some annotated 'Truth' (genes_anno) on a specified region.
% The metric is implemented according to EGASP: (Guigo et al. Genome
% Biology 2006, 7(Suppl I): http://genomebiology.com/content/pdf/gb-2006-7-s1-s2.pdf)
% Additionally we report sensitivity and specificity for tis and stop
% recognition and also for promoter (tss) and cleavage site recognition. 
% for tss and cleave predictions are counted corrected if they are within
% a window of 50,250,1000 nucleotides.  
%
% if genes_anno.mask_on is set to 1, for a given annotated gene, the
% region of this gene is excluded form the evalution. Can be used to look
% on confirmed, unconfirmed orpartially confirmed genes, separately.

verb =1;


if nargin==1
	genes_anno 	= regions.genes_anno; 
	genes_pred 	= regions.genes_pred; 
	genome_info = regions.genome_info; 
	run_local 	= regions.run_local; 
	fn_anno 	= regions.fn_anno; 
	fn_pred 	= regions.fn_pred; 
	fn_output 	= regions.fn_output; 
	regions 	= regions.regions;
	paths
end

if nargin<5
  run_local=1;
end
if 0% ~isequal(unique([genes_pred.chr_num]), unique([genes_anno.chr_num]))
	warning('filtering out annotated genes from contigs where I have no predictions for\n')
	pred_contigs = unique([genes_pred.chr_num]);
	gg = [];
	for j = pred_contigs
		gg = [gg genes_anno([genes_anno.chr_num]==j)];
	end
	genes_anno = gg;
	clear gg pred_contigs
end
if 0%~isequal(unique([genes_pred.chr_num]), unique([genes_anno.chr_num]))
	warning('filtering out predicted genes from contigs where I have no annotations for\n')
	anno_contigs = unique([genes_anno.chr_num]);
	gg = [];
	for j = anno_contigs
		gg = [gg genes_pred([genes_pred.chr_num]==j)];
	end
	genes_pred = gg;
	clear gg anno_contigs 
end

if length(unique([genes_anno.id]))~=length(genes_anno)
	warning('anno ids not unique, assigning new ids');
	for j = 1:length(genes_anno)
		genes_anno(j).id = j;
	end
end

%for cleave and tss
end_tol = [10 20 50 250 1000];
filter_utr_anno=10;
filter_utr_pred=0;
[res, signals] = generate_res;
anno_str = [genes_anno.strand];
pred_str = [genes_pred.strand];
						 
anno_chr = [genes_anno.chr_num];
pred_chr = [genes_pred.chr_num];
num_anno_covered = 0;
num_pred_covered = 0;
id_anno_covered = [];
id_pred_covered = [];

if ~isfield(genes_anno, 'tis'),
  genes_anno(end).tis = [] ;
end ;
if ~isfield(genes_anno, 'cdsStop'),
  genes_anno(end).cdsStop = [] ;
end ;
if ~isfield(genes_anno, 'tss'),
  genes_anno(end).tss = [] ;
end ;
if ~isfield(genes_anno, 'cleave'),
  genes_anno(end).cleave = [] ;
end ;

% evaluation of predictions
Map = [];
fprintf('start evaluatiuon of prediction...\n')
jobinfo = rproc_empty(0) ;
job_num = 0;
done = 0;

%options.addpaths = {'~/svn/projects/genefinding/eval_lsl/','~/svn/projects/genefinding/'} ;
MEMREQ = 4000;
time_req = 60*40;

%regions = regions([regions.chr_num]==12&[regions.strand]=='-');
%regions = regions([regions.strand]=='+');
%regions = regions(9:10);

num_res=0;
while ~done
  for r_idx = 1:length(regions),
    if isfield(regions(r_idx),'strand')
      strands=regions(r_idx).strand;
    else
      strands='+-';
    end
    for s=strands
      fprintf('region %i / %i, strand: %s \n',r_idx, length(regions), s)
      region = regions(r_idx);
      if isfield(region, 'chr_num'),
        chr = region.chr_num;%strmatch(region.chr,genome_info.contig_names,'exact');
      else
        chr = strmatch(region.chr,genome_info.contig_names,'exact');
      end ;

      idx = find(anno_str==s & anno_chr==chr & [genes_anno.start]<=region.stop &[genes_anno.stop]>=region.start);
      anno = genes_anno(idx);
      idx = find(pred_str==s & pred_chr==chr & [genes_pred.start]<=region.stop &[genes_pred.stop]>=region.start);
      pred = genes_pred(idx);    
      
      %if isempty(pred), 
      %  warning('empty prediction') ; 
      %  [tmp,idx1,idx2]=intersect([anno.id], [genes_anno.id]) ;
      %  genes_anno(idx2) = [] ;
      %  anno_str = [genes_anno.strand];
      %  anno_chr = [genes_anno.chr_num];
      %  continue ; 
      %end ;
      
      fprintf('number of annotated genes in region %i%s:\t %i\n',r_idx,s,length(anno));
      fprintf('number of predicted genes in region %i%s:\t %i\n',r_idx,s,length(pred));
      mask_regions =[];
      if isfield(anno,'mask_on')
        idx = find([anno.mask_on]);
        num=0;
        for j=idx
          num=num+1;
          mask_regions(:,num) = [anno(j).start  anno(j).stop]; 
        end
      end
      % keyboard
      if run_local
        fprintf('start mapping annotation...\n')
        %%% annotation
        [A.nuc,A.exons,A.transcripts,A.utr5, A.utr3,A.signals] = map_genes2region(anno,region,mask_regions,filter_utr_anno);
        
        %%% prediction
        fprintf('start mapping predictions...\n')
        [P.nuc,P.exons,P.transcripts,P.utr5, P.utr3,P.signals] = map_genes2region(pred,region,mask_regions,filter_utr_pred);
        done=1;
      else
        clear A
        A.region = region;
        A.mask_regions = mask_regions;
        A.fn = sprintf('%s_%i%s.mat',fn_anno,r_idx,s);
        
        if ~fexist(A.fn)
          options.identifier = ['Map-a' num2str(r_idx),s] ;
          job_num = job_num+1;
          A.genes = anno;
          A.filter_utr = filter_utr_anno;
          jobinfo(job_num) = rproc('map_genes2region_rproc',A,MEMREQ,options,time_req);
          done=0;
        end
        clear P
        P.fn = sprintf('%s_%i%s.mat',fn_pred,r_idx,s);
        P.region = region;
        P.mask_regions = mask_regions;
        if ~fexist(P.fn)
          options.identifier = ['Map-p' num2str(r_idx),s] ;
          job_num = job_num+1;
          P.genes = pred;
          P.filter_utr = filter_utr_pred;
          % map_genes2region_rproc(P)
          jobinfo(job_num) = rproc('map_genes2region_rproc',P,MEMREQ,options,time_req);
          done=0;
        end
        if job_num == 0
          done=1;
        end
        clear A P
      end %if run_local
      % keyboard
      if  done
        if ~run_local
          % [A.nuc,A.exons,A.transcripts,A.utr5, A.utr3,A.signals] = map_genes2region(anno,region,mask_regions,filter_utr_anno);
          % fprintf('loading...\n')
          fn = sprintf('%s_%i%s.mat',fn_anno,r_idx,s);
          A = load(fn);
          % assert(isequal(A,AA))
          % [P.nuc,P.exons,P.transcripts,P.utr5, P.utr3,P.signals] = map_genes2region(pred,region,mask_regions,filter_utr_pred);        
          fn = sprintf('%s_%i%s.mat',fn_pred,r_idx,s);
          P = load(fn);
          % assert(isequal(P,PP))
        end
        % fprintf('count...nucleotides\n')
        %%count observed, predicted, correctly predicted segments
        [res.nucleotides, res.cds_nucleotides] = cnt_obsv_pred_corr(res.nucleotides, res.cds_nucleotides,A.nuc,P.nuc);
        % fprintf('count...exons\n')
        [res.exons, res.cds_exons]             = cnt_obsv_pred_corr(res.exons, res.cds_exons,A.exons,P.exons);
        % fprintf('count...transcripts\n')
        [res.transcripts, res.cds_transcripts] = cnt_obsv_pred_corr(res.transcripts, res.cds_transcripts,A.transcripts,P.transcripts);
        % fprintf('count...signals\n')
        signals  = cnt_obsv_pred_corr(signals, [],A.signals,P.signals,s,end_tol);
        % fprintf('count...genes\n')
        %map genes together
        [res.genes,map] = compare_genes(res.genes,anno,pred);
        
        num_res=num_res+1;
        [RES, SIGNALS] = generate_res;
         % fprintf('count...nucleotides\n')
        %%count observed, predicted, correctly predicted segments
        [RES.nucleotides, RES.cds_nucleotides] = cnt_obsv_pred_corr(RES.nucleotides, RES.cds_nucleotides,A.nuc,P.nuc);
        % fprintf('count...exons\n')
        [RES.exons, RES.cds_exons]             = cnt_obsv_pred_corr(RES.exons, RES.cds_exons,A.exons,P.exons);
        % fprintf('count...transcripts\n')
        [RES.transcripts, RES.cds_transcripts] = cnt_obsv_pred_corr(RES.transcripts, RES.cds_transcripts,A.transcripts,P.transcripts);
        [RES.genes,map] = compare_genes(RES.genes,anno,pred);
        R(num_res) = compute_SN_SP(RES,signals,end_tol);
      else
        R=[];
        % Map(:,end+[1:size(map,2)])=map;
      end 
      
    end %loop over strand
  end % loop over r_idx
  % [jobinfo,num_crashed] = rproc_wait(jobinfo) ;
  done=1
end %while not done

% keyboard
fprintf('\n loop over regions done\n')


% Compute sensitivity and specificity
res = compute_SN_SP(res,signals,end_tol);



% keyboard

if exist('fn_output','var') 
  fd = fopen(fn_output, 'w+') ;
  if fd<1 && verb
	fds = 1;
  elseif verb
	fds = [fd 1];
  end
else
  fds = 1;
end
% Print evaluation summary
for fd = fds
  fprintf(fd,'Evaluation\n');
  fprintf(fd,'  on coding nucleotide level:\n');
  fprintf(fd,'    sensitivity: %2.2f%%\n', 100*res.cds_nucleotides.SN);
  fprintf(fd,'    specificity: %2.2f%%\n\n', 100*res.cds_nucleotides.SP);
  
  fprintf(fd,'  on coding exon level:\n');
  fprintf(fd,'    sensitivity: %2.2f%% (missing exons: %2.2f%%)\n', ...
	  100*res.cds_exons.SN, 100*res.cds_exons.ME/res.cds_exons.num_obsv);
  fprintf(fd,'    specificity: %2.2f%% (wrong exons:   %2.2f%%)\n\n', ...
	  100*res.cds_exons.SP, 100*res.cds_exons.WE/res.cds_exons.num_pred);
  
  fprintf(fd,'  on coding transcript level(F:%2.2f):\n',100*fscore(res.cds_transcripts.SN,res.cds_transcripts.SP));
  fprintf(fd,'    sensitivity: %2.2f%% \n',100*res.cds_transcripts.SN);
  fprintf(fd,'    specificity: %2.2f%% \n\n' , 100*res.cds_transcripts.SP);
  
  fprintf(fd,'  on gene level:\n');
  fprintf(fd,'    sensitivity: %2.2f%% \n',100*res.genes.SN);
  fprintf(fd,'    specificity: %2.2f%% \n\n', 100*res.genes.SP);
  
  fprintf(fd,'    missing genes: %i (%2.2f%%)\n',res.genes.missing,100*res.genes.missing/res.genes.num_obsv);
  fprintf(fd,'    wrong genes: %i (%2.2f%%)\n\n',res.genes.wrong,100*res.genes.wrong/res.genes.num_pred);
  
  fprintf(fd,'  TIS:\n');
  fprintf(fd,'    sensitivity: %2.2f%%\n', 100*res.tis.SN);
  fprintf(fd,'    specificity: %2.2f%%\n\n', 100*res.tis.SP);
  
  fprintf(fd,'  STOP:\n');
  fprintf(fd,'    sensitivity: %2.2f%%\n', 100*res.cdsStop.SN);
  fprintf(fd,'    specificity: %2.2f%%\n\n', 100*res.cdsStop.SP);
 
  if fd>1
  	fclose(fd);
  end
end

return
