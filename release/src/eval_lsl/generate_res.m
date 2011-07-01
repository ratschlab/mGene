function [res signals] = generate_res


res.nucleotides.num_obsv = 0;
res.nucleotides.num_pred = 0;
res.nucleotides.num_corr = 0;

res.cds_nucleotides.num_obsv = 0;
res.cds_nucleotides.num_pred = 0;
res.cds_nucleotides.num_corr = 0;

res.exons.num_obsv = 0;
res.exons.num_pred = 0;
res.exons.num_corr = 0;
res.exons.ME = 0; % missing exons (no overlap to pred)
res.exons.WE = 0; % wrong exons (no overlap to anno)
res.exons.ME_gene_id =[];
res.exons.WE_gene_id =[];

res.cds_exons.num_obsv = 0;
res.cds_exons.num_pred = 0;
res.cds_exons.num_corr = 0;
res.cds_exons.ME = 0; % missing exons (no overlap to pred)
res.cds_exons.WE = 0; % wrong exons (no overlap to anno)
res.cds_exons.ME_gene_id =[];
res.cds_exons.WE_gene_id =[];

res.transcripts.num_obsv = 0;
res.transcripts.num_pred = 0;
res.transcripts.num_corr = 0;

res.cds_transcripts.num_obsv = 0;
res.cds_transcripts.num_pred = 0;
res.cds_transcripts.num_corr = 0;

res.genes.num_obsv = 0;
res.genes.num_pred = 0;
res.genes.num_corr = 0;
res.genes.wrong    = 0;
res.genes.wrong_ids= [];
res.genes.missing  = 0;
res.genes.missing_ids = [];
res.genes.missing_conf = 0;
res.genes.missing_partconf = 0;
res.genes.missing_unconf = 0;
res.genes.cut = 0;
res.genes.merged = 0;


status=inf;

%%%%%
signals.acc.num_obsv = 0;
signals.acc.num_pred = 0;
signals.acc.num_corr = 0;

signals.don.num_obsv = 0;
signals.don.num_pred = 0;
signals.don.num_corr = 0;

signals.tis.num_obsv = 0;
signals.tis.num_pred = 0;
signals.tis.num_corr = 0;
signals.tis.dist = [];

signals.cdsStop.num_obsv = 0;
signals.cdsStop.num_pred = 0;
signals.cdsStop.num_corr = 0;
signals.cdsStop.dist = [];

signals.tss.num_obsv = 0;
signals.tss.num_pred = 0;
signals.tss.num_corr = 0;

signals.tss.num_corr1000 = 0;
signals.tss.num_corr250 = 0;
signals.tss.num_corr50 = 0;
signals.tss.num_corr20 = 0;
signals.tss.num_corr10 = 0;

signals.tss.dist = [];
signals.tss.utr5p_length_obsv = [];
signals.tss.utr5p_length_pred = [];

signals.cleave.num_obsv = 0;
signals.cleave.num_pred = 0;
signals.cleave.num_corr = 0;
signals.cleave.num_corr1000 = 0;
signals.cleave.num_corr250 = 0;
signals.cleave.num_corr50 = 0;
signals.cleave.num_corr20 = 0;
signals.cleave.num_corr10 = 0;

signals.cleave.dist = [];
signals.cleave.utr3p_length_obsv = [];
signals.cleave.utr3p_length_pred = [];
