function Content=content_template(name)

Content.name = name;
Content.consensus = [];
Content.lwin_big = 0;
Content.rwin_big = 0;
%Content.order = [3 4 5];
Content.label_fct = 'get_cand_contents';
Content.load_fct = 'load_content_data';

Content.filter_label.train.USE_ALL = 1;
Content.filter_label.test.USE_ALL = 1;
Content.filter_label.eval.USE_ALL = 1;
Content.Conf_names = {} ;
Content.C = [1]; 

Content.wordlen = { };
Content.stepping = 1;
Content.offset = 0;

Content.export_settings.conf_cum_thresh = 0.5 ;
Content.export_settings.resolution = 1 ;


%PAR.Contents.type_descrs = { 'intergenic', 'utr_exon', 'cds_exon','intron', 'intercistronic' };
%PAR.Contents.type_abbrs = { 'ige' 'utr' 'cex' 'ino' 'ics' };
%PAR.Contents.num_types = length(PAR.Contents.type_descrs); 


