function [transcript_hash_names, transcript_hash_idx]= init_hash_search(transcript_names, new_start, transcript_hash_names, transcript_hash_idx) ;
% [transcript_hash_names, transcript_hash_idx]= init_hash_search(transcript_names, new_start, transcript_hash_names, transcript_hash_idx) ;

if nargin<=1,
  transcript_hash_names = {} ;
  transcript_hash_names{65536}={} ;
  transcript_hash_idx = {} ;
  transcript_hash_idx{65536} = [] ;
  new_start = 1 ;
end ;
for i=new_start:length(transcript_names),
        if mod(i,1000)==0, fprintf('generated %i hash keys\n', i) ; end ;
	checksum=crc16(transcript_names{i})+1 ;
        transcript_hash_names{checksum}{end+1} = transcript_names{i} ;
        transcript_hash_idx{checksum}(end+1) = i ;
end ;
