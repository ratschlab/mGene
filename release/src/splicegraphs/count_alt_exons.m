function count_alt_exons(data)

addpath detect_altsplice/
addpath ~/svn/tools/utils/


count_exons = 1;
count_genes = 0;
count_nematode = 0;

organismlistfile = data.organismlistfile;
outfile = data.outfile;
summaryfile = data.summaryfile;
sourcefile = data.sourcefile;
fd_rate = data.fd_rate;


orglist = textread(organismlistfile,'%s')
numorgs = length(orglist);


frac_exon = zeros(numorgs,1);
frac_intron = zeros(numorgs,1);
frac_5prime = zeros(numorgs,1);
frac_3prime = zeros(numorgs,1);

%basedir = '/local/cong/Data/altsplice';
basedir = '/fml/ag-raetsch/share/projects/altsplicedata';

fid = fopen(outfile,'w');

if count_exons
  fprintf(1,'\n\nCounting Exons\n\n');
  fprintf(fid,'\n\nCounting Exons\n\n');

  for ix = 1:numorgs
    fprintf(1,'\n\n%%%% %s %%%%\n\n',orglist{ix});
    fprintf(fid,'\n\n%%%% %s %%%%\n\n',orglist{ix});
    if strcmp(sourcefile,'confirmed')
      genesfile = sprintf('%s/%s/confirmed_sequences.mat',basedir,orglist{ix});
    elseif strcmp(sourcefile,'pred')
      genesfile = sprintf('%s/%s/prediction/cand_splicegraph_%02d.mat',basedir,orglist{ix},fd_rate);
    else
      fprintf('sourcefile must be "confirmed" or "pred"\n');
    end
    genes = load_cell(genesfile,'genes');
    [frac_exon(ix), frac_intron(ix), frac_5prime(ix), frac_3prime(ix),nums{ix},totals{ix}] = exon_counter(genes,fid);
  end
  
  fprintf(1,'\n\n\n');
  fprintf(fid,'\n\n\n');

  for ix = 1:numorgs
    fprintf(fid,'%s,%f,%f,%f,%f\n',orglist{ix},frac_exon(ix),frac_intron(ix),frac_5prime(ix),frac_3prime(ix));
    fprintf(1,'%s,%f,%f,%f,%f\n',orglist{ix},frac_exon(ix),frac_intron(ix),frac_5prime(ix),frac_3prime(ix));
  end


  fprintf(fid,'\n\n');
  fprintf(1,'\n\n');

  fid2 = fopen(summaryfile,'w');
  for ix = 1:numorgs
    fprintf(fid2,'%s,%d,%d,%d,%d\n',orglist{ix},nums{ix}.exon,nums{ix}.intron,nums{ix}.prime5,nums{ix}.prime3);
    fprintf(fid,'%s,%d,%d,%d,%d\n',orglist{ix},nums{ix}.exon,nums{ix}.intron,nums{ix}.prime5,nums{ix}.prime3);
    fprintf(1,'%s,%d,%d,%d,%d\n',orglist{ix},nums{ix}.exon,nums{ix}.intron,nums{ix}.prime5,nums{ix}.prime3);
    fprintf(fid2,'%s,%d,%d,%d,%d\n',orglist{ix},totals{ix}.exon,totals{ix}.intron,totals{ix}.prime5,totals{ix}.prime3);
    fprintf(fid,'%s,%d,%d,%d,%d\n',orglist{ix},totals{ix}.exon,totals{ix}.intron,totals{ix}.prime5,totals{ix}.prime3);
    fprintf(1,'%s,%d,%d,%d,%d\n',orglist{ix},totals{ix}.exon,totals{ix}.intron,totals{ix}.prime5,totals{ix}.prime3);
  end
  fclose(fid2);
end


  



if count_genes
  fprintf(1,'\n\nCounting Genes\n\n');
  fprintf(fid,'\n\nCounting Genes\n\n');

  for ix = 1:numorgs
    fprintf(1,'\n\n%%%% %s %%%%\n\n',orglist{ix});
    fprintf(fid,'\n\n%%%% %s %%%%\n\n',orglist{ix});
    if strcmp(sourcefile,'confirmed')
      genesfile = sprintf('%s/%s/confirmed_sequences.mat',basedir,orglist{ix});
    elseif strcmp(sourcefile,'pred')
      genesfile = sprintf('%s/%s/prediction/cand_splicegraph_%02d.mat',basedir,orglist{ix},fd_rate);
    else
      fprintf('sourcefile must be "confirmed" or "pred"\n');
    end
    load(genesfile,'genes');
    [frac_exon(ix),frac_intron(ix),frac_5prime(ix),frac_3prime(ix),nums{ix},totals{ix}] = genes_counter(genes,fid);
  end

  fprintf(1,'\n\n\n');
  fprintf(fid,'\n\n\n');
  
  for ix = 1:numorgs
    fprintf(fid,'%s,%f,%f,%f,%f\n',orglist{ix},frac_exon(ix),frac_intron(ix),frac_5prime(ix),frac_3prime(ix));
    fprintf(1,'%s,%f,%f,%f,%f\n',orglist{ix},frac_exon(ix),frac_intron(ix),frac_5prime(ix),frac_3prime(ix));
  end


  fprintf(fid,'\n\n');
  fprintf(1,'\n\n');

  fid2 = fopen(summaryfile,'w');
  for ix = 1:numorgs
    fprintf(fid,'%s,%d,%d,%d,%d\n',orglist{ix},nums{ix}.exon,nums{ix}.intron,nums{ix}.prime5,nums{ix}.prime3);
    fprintf(1,'%s,%d,%d,%d,%d\n',orglist{ix},nums{ix}.exon,nums{ix}.intron,nums{ix}.prime5,nums{ix}.prime3);
    fprintf(fid,'%s,%d,%d,%d,%d\n',orglist{ix},totals{ix}.exon,totals{ix}.intron,totals{ix}.prime5,totals{ix}.prime3);
    fprintf(1,'%s,%d,%d,%d,%d\n',orglist{ix},totals{ix}.exon,totals{ix}.intron,totals{ix}.prime5,totals{ix}.prime3);
  end
  fclose(fid2);
end



fclose(fid);






