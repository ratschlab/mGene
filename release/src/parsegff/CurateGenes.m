function [genes] = CurateGenes(genome_info_path, gene_models_file)
    genome_info_file = strcat(genome_info_path, '/genome.config');

    disp('----------------------------------');
    disp('Step 1: Reading Genome Info in GIO');
    disp('----------------------------------'); 
    fprintf('\nReading Genome Info in GIO : %s\n', genome_info_path);
    genome_info = init_genome(genome_info_file);
    fprintf('Done.\n\n');
    
    disp('-------------------------------------');
    disp('Step 2: Reading Gene models from GFF3');
    disp('-------------------------------------');
    fprintf('\nReading Gene models from : %s\n', gene_models_file);
    load(gene_models_file);
    fprintf('\nRead %i Gene models\n', length(genes));
    fprintf('Done.\n\n');
    
    disp('-----------------------------------------------------');
    disp('Step 3: Validating Gene models and Genome Info in GIO');
    disp('-----------------------------------------------------');
    
    gff_chr_names = unique({genes.chr});
    unmap_chr = setdiff(gff_chr_names, genome_info.contig_names);
    if ~isempty(unmap_chr),
        error('Warning: The Gene models from the following chromosomes/contigs are not present in Genome Info GIO:\n'); 
        for i=1:length(unmap_chr),
            fprintf(1, '     %s\n', unmap_chr{i});
        end;
        quit; 
    end;
    fprintf('Done.\n\n');
    
    % adding chr_num to the gene models
    for i=1:length(genes)
        genes(i).chr_num = strmatch(upper(genes(i).chr), upper(genome_info.contig_names), 'exact');
    end; 
    
    disp('-----------------------------------------------------------');
    disp('Step 4: Checking Translation Initiation Site in Gene models');
    disp('-----------------------------------------------------------');
    genes = correct_tis_stop(genes, genome_info);
    fprintf('Done.\n\n');
    % checking parsed gene models
    which_checks.exons_sorted = 1;
    which_checks.intron_length = 1;
    which_checks.splicesites = 1;
    which_checks.orf = 1;
    which_checks.gene_length = 1;
    which_checks.graph = 0;
    which_checks.transacc = 0;
    which_checks.complete = 1;
    disp('----------------------------');
    disp('Step 5: Checking Gene models');
    disp('----------------------------');
    genes = check_genes(genes, genome_info, which_checks);
    disp('-----------------------------');
    disp('Step 6: Filtering Gene models');
    disp('-----------------------------');
    genes = filter_invalid_genes(genes, genome_info);

    % curating gene models
    disp('----------------------------');
    disp('Step 7: Curating Gene models');
    disp('----------------------------');
    genes = curate_anno(genes, genome_info);
 
    % checking gene models
    genes = check_genes(genes, genome_info, which_checks);
    
    % build splice-graph
    disp('----------------------------');
    disp('Step 8: Building splicegraph');
    disp('----------------------------');
    genes = update_splicegraph(genes);    
    fprintf('Done.\n\n');
    disp('---------------------------------------------------');
    disp('Step 9: Generating Annotation Gene Structure object');
    disp('---------------------------------------------------');

    save('-v7', gene_models_file, 'genes');

    fprintf('\nAnnotation with %i genes saved at %s\n', length(genes), gene_models_file);
    fprintf('Done.\n');
