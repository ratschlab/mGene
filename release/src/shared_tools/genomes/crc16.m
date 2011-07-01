function checksum = crc16(String)  
% checksum = crc16(String)  
 
checksum = sg('crc', String) ;
checksum = uint32(checksum) ;
checksum = bitand(checksum, 65535) + floor(checksum/65535) ;
checksum = bitand(checksum, 65535) ;

