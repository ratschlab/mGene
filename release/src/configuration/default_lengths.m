function lengths_range = default_lengths

lengths_range.intergenic = [0, 21000] ;
%lengths_range.intergenictrans = [0, 21000] ;
lengths_range.utr5exon   = [0, 2000] ;
% lengths_range.transexon  = [1, 1000] ;
lengths_range.single_cds_exon = [100, 8000] ; %%%100, from don => acc+stop
lengths_range.first_cds_exon  = [3, 8000] ;
lengths_range.middle_cds_exon   = [3, 8000] ;
lengths_range.last_cds_exon   = [1, 8000] ;
lengths_range.utr3exon    = [1, 2000] ;
lengths_range.polya_tail  = [0, 250] ;
% lengths_range.intercistronic = [1, 5000] ;
lengths_range.intergenic_long = [10500-50, 10500+50] ; % due to irregular spacing 
%lengths_range.intron      = [30, 100000] ;%%%30, from don => acc+stop
lengths_range.intron      = [30, 20000] ;%%%30, from don => acc+stop
lengths_range.rna_seq_polya = [0, 2000];
lengths_range.first_nc_exon = [0, 8000];
lengths_range.middle_nc_exon = [0, 8000];
lengths_range.last_nc_exon = [0, 8000];
