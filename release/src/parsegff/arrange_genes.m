function genes = arrange_genes(filename)
    try
        l = load(filename);
		if isfield(l, 'arrange_genes_done')
			fprintf('arrange_genes was already applied to this file: %s', filename);
			genes = l.genes;
			return
		end
    catch 
		fprintf('could not load pre-processed annotation in .mat format %s\n', filename)
        quit; 
    end
    x = 1;
    for i=1:length(l.gene_models) 
        gene = l.gene_models{x};  
        gene.alias = '';
        % arranging transcript_id
        gene.transcripts = reshape(gene.transcripts, 1, length(gene.transcripts));
        % arranging transcript_status
        gene.transcript_status = double(gene.transcript_status);
        gene.transcript_status = reshape(gene.transcript_status, 1, length(gene.transcript_status));
        % reshaping exons
        for j=1:length(gene.exons)
            gene.exons{j} = double(gene.exons{j});
        end    
        gene.exons = reshape(gene.exons, 1, length(gene.exons));
        % exons_confirmed
        ex = {};
        for ec=1:length(gene.transcripts)
            ex{1,ec} = [];
        end
        gene.exons_confirmed = ex;
        % reshaping cds_exons
        for j=1:length(gene.cds_exons)
            gene.cds_exons{j} = double(gene.cds_exons{j});
        end
        gene.cds_exons = reshape(gene.cds_exons, 1, length(gene.cds_exons));
        % reshaping utr5_exons
        for j=1:length(gene.utr5_exons)
            gene.utr5_exons{j} = double(gene.utr5_exons{j}); 
        end 
        gene.utr5_exons = reshape(gene.utr5_exons, 1, length(gene.utr5_exons));
        % reshaping utr3_exons
        for j=1:length(gene.utr3_exons)
            gene.utr3_exons{j} = double(gene.utr3_exons{j});
        end
        gene.utr3_exons = reshape(gene.utr3_exons, 1, length(gene.utr3_exons));
        % reshaping tis
        for j=1:length(gene.tis)
            gene.tis{j} = double(gene.tis{j});
        end
        gene.tis = reshape(gene.tis, 1, length(gene.tis));
        % reshaping cdsStop
        for j=1:length(gene.cdsStop)
            gene.cdsStop{j} = double(gene.cdsStop{j});
        end
        gene.cdsStop = reshape(gene.cdsStop, 1, length(gene.cdsStop));
        for j=1:length(gene.tss)
            gene.tss{j} = double(gene.tss{j});
        end
        gene.tss = reshape(gene.tss, 1, length(gene.tss));
        for j=1:length(gene.cleave)
            gene.cleave{j} = double(gene.cleave{j});
        end
        gene.cleave = reshape(gene.cleave, 1, length(gene.cleave));
        gene.id = double(gene.id);
        gene.start = double(gene.start);
        gene.stop = double(gene.stop);
		genes(length(l.gene_models)-i+1) = gene; 
        x = x + 1;
    end
	arrange_genes_done=1;
    save('-v7', filename, 'genes', 'arrange_genes_done');
    clear l
