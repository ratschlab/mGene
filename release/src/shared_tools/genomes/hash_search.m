function idx = hash_search(transcript_hash_names, transcript_hash_idx, transcript_name) ;
% idx = hash_search(transcript_hash_names, transcript_hash_idx, transcript_name) ;

checksum=crc16(transcript_name)+1 ;
idx_ = strmatch(transcript_name, transcript_hash_names{checksum}, 'exact') ;
idx = transcript_hash_idx{checksum}(idx_) ;
