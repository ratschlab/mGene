% copy of statistics(info, origpath, directory)
% except with the setup information removed.
% assume genes structure exists.





%addpath(sprintf('%s/utils', origpath)) ;
%addpath(sprintf('%s/splicegraphs', origpath)) ;
%addpath(sprintf('%s/splicegraphs/detect_altsplice', origpath)) ;
%addpath(sprintf('%s/gff', origpath)) ;
%addpath(sprintf('%s/gff/read_sequences', origpath)) ;
%addpath /fml/ag-raetsch/share/software/matlab_tools/utils
%addpath /fml/ag-raetsch/share/software/matlab_tools/rproc/

%directory='statistics' ;

%%%%%%%%%%%%%%%%%%%%%
%%% Genome.Config %%%
%%%%%%%%%%%%%%%%%%%%%
%genome_info = init_genome(info) ;
%fprintf('directory: %s\n', genome_info.basedir) ;
%%%%genome_info = init_genome('/fml/ag-raetsch/share/projects/splicing/test3/genome.config');


%%%%%%%%%%%%%%%
%%% Loading %%%
%%%%%%%%%%%%%%%
%load(sprintf('%s/confirmed_sequences.mat', genome_info.basedir), 'genes') ;
%fprintf('loaded confirmed_sequences.mat\n\n') ;


%%%%%%%%%%%%%%%%
%%% Saved in %%%
%%%%%%%%%%%%%%%%
%fprintf('results in: %s/%s\n', genome_info.basedir, directory) ;
%file_name='alternative_numbers.txt' ;
%fid=fopen(sprintf('%s/%s/%s', genome_info.basedir, directory, file_name),'w+') ;
%fprintf(fid, 'directory: %s\n', genome_info.basedir) ;
%fprintf(fid, 'loaded confirmed_sequences.mat') ;


%%%%%%%%%%%%%%%%%
%%% ALT-CONST %%%
%%%%%%%%%%%%%%%%%

if ~exist('fid'), fid=1 ; end ;

genes=alt_const(genes) ;
fprintf(fid,'\n\nTotal genes with alternative isoforms:\t\t\t\t%d\n',...
        sum([genes(:).is_alt]));
fprintf(fid,'Total genes alternatively spliced:\t\t\t\t%d\n',...
        sum([genes(:).is_alt_spliced]));
fprintf(fid,'Total constitutively spliced:\t\t\t\t\t%d\n',...
        sum(~[genes(:).is_alt_spliced]));



[idx_alt_tstart, exon_alt_tstart] = detect_alttstart(genes) ;
fprintf(fid,'\nTotal alternative transcription starts:\t\t\t\t%d\n',...
        length(idx_alt_tstart));


[idx_alt_tend, exon_alt_tend] = detect_alttend(genes) ;
fprintf(fid,'\nTotal alternative transcription ends:\t\t\t\t%d\n',...
        length(idx_alt_tend));

idx_alt=[] ;
idx_test=[] ;
for i=1:length(genes)
  if genes(i).is_alt_spliced
    idx_alt=[idx_alt, i] ;
  end
  if genes(i).is_alt
    idx_test=[idx_test, i] ;
  end
end

genes_test=genes(idx_test) ;
[idx_alt_tstart_test, exon_alt_tstart_test] = detect_alttstart(genes_test) ;
[idx_alt_tend_test, exon_alt_tend_test] = detect_alttend(genes_test) ;
assert(length(idx_alt_tstart) == length(idx_alt_tstart_test)) ;
assert(length(idx_alt_tend) == length(idx_alt_tend_test)) ;


%%%%%%%%%%%%%%%%%%
%%% Statistics %%%
%%%%%%%%%%%%%%%%%%


%%% count new exons created by infer %%%
[genes,new_exon_count]=detect_newexons(genes, idx_alt) ;
fprintf(fid,'\nNumber of new exons:\t\t\t\t\t\t%d\n', ...
        new_exon_count);



% detect exon skips %
[idx_exon_skips, exon_exon_skips] = detect_exonskips(genes, idx_alt);
fprintf(fid,'\n\nNumber of single exon skips:\t\t\t\t\t%d\n', ...
        length(idx_exon_skips));
  


%%% detect alternative 5 and 3 prime sites %%%
[idx_alt_5prime,exon_alt_5prime, idx_alt_3prime,exon_alt_3prime] = ...
    detect_altprime(genes, idx_alt);
fprintf(fid,'\n\nNumber of alternative 5 prime sites:\t\t\t\t%d\n',...
	length(idx_alt_5prime));
fprintf(fid,'Number of alternative 3 prime sites:\t\t\t\t%d\n',...
	length(idx_alt_3prime));



%%% detect intron retentions %%%
[idx_intron_reten,intron_intron_reten] = detect_intronreten(genes, idx_alt) ;
fprintf(fid,'\n\nNumber of intron retentions:\t\t\t\t\t%d\n', ...
	length(idx_intron_reten));



%%% detect XOR exons %%%
[idx_xor_exons, exon_xor_exons] = detect_xorexons(genes, idx_alt) ;
fprintf(fid,'\n\nNumber of XOR exons:\t\t\t\t\t\t%d\n',...
	length(idx_xor_exons));



idx_multiple_skips=[] ;
%%% detect multiple exon skips %%%
if 1
  [idx_multiple_skips, exon_multiple_skips] = ...
      detect_multipleskips(genes, idx_alt) ;
  fprintf(fid,'\nNumber of multiple exon skips:\t\t\t\t\t%d\n',...
	length(idx_multiple_skips));
end



%%% detect incomplete exons in intronic regions %%%
[idx_exon_intron5, exon_exon_intron5, idx_exon_intron3, exon_exon_intron3] = ...
    detect_exonintron(genes, idx_alt) ; 
fprintf(fid,'\n\nNumber of incomplete exons in intronic regions (5ps):\t\t%d\n',...
	length(idx_exon_intron5));
fprintf(fid,'Number of incomplete exons in intronic regions (3ps):\t\t%d\n',...
	length(idx_exon_intron3));




%%% detect alternative introns %%%
[idx_alt_intron, introns_alt_intron] = detect_altintrons(genes, [1:length(genes)]) ;
fprintf(fid,'\n\nNumber of alternative introns:\t\t\t\t\t%d\n',...
	length(idx_alt_intron));



%%% detect undetermined splicing events %%%            
idx_unknown = [];
for ix=idx_alt
  if (mod(ix,50)==0)
    fprintf(1,'.');
  end
  if (isempty(find(ix==idx_exon_skips))) && ...
        (isempty(find(ix==idx_alt_5prime))) && ...
        (isempty(find(ix==idx_alt_3prime))) && ...
        (isempty(find(ix==idx_intron_reten))) && ...
        (isempty(find(ix==idx_alt_intron))) && ...  
        (isempty(find(ix==idx_xor_exons))) && ...
        (isempty(find(ix==idx_multiple_skips))) && ...
        (isempty(find(ix==idx_exon_intron5))) && ...
        (isempty(find(ix==idx_exon_intron3))),
    idx_unknown = [idx_unknown,ix];
  end
end
fprintf(1,'\n\nNumber of undetermined splicing events:\t\t\t\t%d\n', length(idx_unknown));
fprintf(fid,'\n\nNumber of undetermined splicing events:\t\t\t\t%d\n', length(idx_unknown));

if fid~=1,
  fclose(fid) ;
end ;

