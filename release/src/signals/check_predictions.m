function all_good = check_predictions(fn_genome_config,fn_pred,Signal_name,resolution)
% check_predictions(fn_pred,fn_genome_config)



P.fn_genome_config = fn_genome_config;
P.fn_pred = fn_pred;
P.resolution = resolution;
P.Signal_name = Signal_name;

genome_info = init_genome(fn_genome_config);

all_good = zeros(2, length(genome_info.contig_names));
for c = 1:length(genome_info.contig_names)
	  num = 0;
   for s = '+-'
     num = num+1;
     filename = sprintf('%scontig_%i%s', fn_pred, c, s);
     if fexist([filename '.pos'])&& fexist([filename '.output'])
       P.chrom = c;
       P.strand = s; 
       all_good(num,c) = check_predictions_helper(P) ;
       fprintf('check_predictions: strand: %s, chrom: %i, ok: %i',s, c,all_good);
     else
       fprintf('file missing: %s',filename)
     end 
   end  
end 
