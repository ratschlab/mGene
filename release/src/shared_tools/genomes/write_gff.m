function write_gff3(genes,gff_fname,source,status)

% write_gff3(genes,gff_fname,source)

  
if nargin<3
  warning('source empty; using source = mGene')
  source = 'mGene';
end
if nargin<4
  warning('status empty; using status = predicted')
  status = 'Predicted';
end

%%% SOFA mapping
% SOFA_short_names = {'gene','mRNA','exon','CDS','intron','polyA_sequence','polyA_site','five_prime_UTR','three_prime_UTR'};
% SOFA_IDs = {'SO:0000704','SO:0000234','SO:0000147','SO:0000316','SO:0000188','SO:0000610','SO:0000553','SO:0000204','SO:0000205'};
% hierachies = [1 2 3 3 3 3 3 3 3];
% status={'Predicted','Partially_confirmed','Confirmed'};
STATS_genes = genes_statistics(genes,source)

% initialize outfile
if fexist(gff_fname)
  fprintf('replacing file %s...\n', gff_fname);
  [fd msg] = fopen(gff_fname, 'w+');
  disp(msg);
else
  fprintf('creating file %s...\n', gff_fname);
  [fd msg] = fopen(gff_fname, 'w+');
  disp(msg);
end

fprintf(fd,'##gff-version 3\n');


cnt_alt = 0;
for g=1:length(genes),
  fprintf('writing gene %i...\r', g);  
  gene = genes(g);
  type  = 'gene';
  score = '.';
  phase = '.';
  attr_str = sprintf('ID=Gene%s' ,gene.name);
  fprintf(fd,'%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n', ...
          gene.chr,source,type,gene.start,gene.stop,score,gene.strand,phase,attr_str);
  for t=1:length(gene.transcripts)
    type  = 'mRNA';
    score = '.';
    phase = '.';
    if gene.strand=='+'
      start = min(gene.exons{t}(:,1));
      stop = max(gene.exons{t}(:,2))-1;
    else
      start = min(gene.exons{t}(:,1))+1;
      stop = max(gene.exons{t}(:,2));
    end

    attr_str = sprintf('ID=Transcript:%s.%i;Parent=gene%s;prediction_status=',gene.name,t,gene.name,status);
    fprintf(fd,'%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n', ...
            gene.chr,source,type,start,stop,score,gene.strand,phase,attr_str);
    exon_type = {'utr5_exons', 'cds_exons', 'utr3_exons'};
    gff_types  = {'five_prime_UTR','CDS','three_prime_UTR'};
    for tt=1:3
      exons=gene.(exon_type{tt}){t};
      for e=1:size(exons,1)
        type  = gff_types{tt};
        score = '.';
        if tt==2
          phase = num2str(exons(e,3));
        else
          phase='.';
        end
        if gene.strand=='+'
          start = exons(e,1);
          stop = exons(e,2)-1;
        else
          start = exons(e,1)+1;
          stop = exons(e,2);
        end
        attr_str = sprintf('ID=%s:%s.%i;Parent=Transcript:%s.%i',gff_types{tt},gene.name,t,gene.name,t);
        fprintf(fd,'%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n', ...
                gene.chr,source,type,start,stop,score,gene.strand,phase,attr_str);
      end %%loop over exons
    end %% loop over types
  end %%% loop over transcripts
end %%% loop over genes


fclose(fd);
[tmp,STATS_gff] = gff3_statistics(gff_fname, source)
types=fieldnames(STATS_genes)
for i=1:length(types)
  assert(STATS_gff.(types{t})==STATS_genes.(types{t}))
end