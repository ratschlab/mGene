function [transcript_hash_names, transcript_hash_idx]= init_hash_search(transcript_names) ;
% [transcript_hash_names, transcript_hash_idx]= init_hash_search(transcript_names) ;

transcript_hash_names = {} ;
transcript_hash_names{65536}={} ;
transcript_hash_idx = {} ;
transcript_hash_idx{65536} = [] ;
for i=1:length(transcript_names),
	checksum=crc16(transcript_names{i}) ;
        transcript_hash_names{checksum}{end+1} = transcript_names{i} ;
        transcript_hash_idx{checksum}(end+1) = i ;
end ;
