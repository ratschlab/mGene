function [organism label_source] = dummy_organism(label_source_orig,signal_names)

[label_source] = organism_template(label_source_orig,signal_names);
  
organism.clade = [];
organism.name = '';
organism.full_name = '';
organism.release = '';
organism.version = '';
organism.genebuilt = '';

