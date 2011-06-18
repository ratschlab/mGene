function lpenv = cplex_license(tmp1, tmp2)

license_file = getenv('ILOG_LICENSE_FILE');
envstr = sprintf('ILOG_LICENSE_FILE=%s', license_file);

global lpenv ;

lpenv = 0 ;
cnt = 0;
while lpenv==0 && cnt<10
	cnt = cnt+1;
	fprintf('\ntrying license file %s\n', strrep(envstr, 'ILOG_LICENSE_FILE=', ''));
	lpenv = cplex_init_quit(0, envstr) ;

	if lpenv==0,
		disp('waiting for cplex license') ;
		pause(60)
	end;
end ;

if lpenv==0
	error('could not get cplex license\n')
end
