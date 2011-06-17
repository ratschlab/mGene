function lpenv = cplex_license(waitflag, num_cpus)

license_file = getenv('ILOG_LICENSE_FILE');
envstr = sprintf('ILOG_LICENSE_FILE=%s', license_file);

global lpenv ;

if nargin<1, waitflag = 0 ; end ;
if nargin<2, num_cpus = 1 ; end ;

lpenv = 0 ;
cnt = 0;
while lpenv==0 && cnt<10
	cnt = cnt+1;
	fprintf('\ntrying license file %s\n', strrep(envstr, 'ILOG_LICENSE_FILE=', ''));
	lpenv = cplex_init_quit(0, envstr) ;

	if lpenv == 0 & ~waitflag,
		break ;
	end ;
	if lpenv==0,
		disp('waiting for cplex license') ;
		pause(60)
	end;
end ;

if lpenv==0
	error('could not get cplex license\n')
end
