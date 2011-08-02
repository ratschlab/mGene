function run_matlab_command(function_name, source_path, varargin)

global SHELL_INTERPRETER_INVOKE; 
SHELL_INTERPRETER_INVOKE=1; 
addpath(source_path); 
paths; 

try, 
	fprintf('calling: \n\t%s(', function_name);
	warning('off')
	for j = 1:length(varargin)
		if j<length(varargin)
			fprintf('''%s'',', varargin{j})
		else
			fprintf('''%s'')\n', varargin{j})
		end
	end
	feval(function_name, varargin{:}); 
catch; 
	err=lasterror;
	fprintf('Error: %s\n', err.message)
	for i=1:length(err.stack) 
		disp(sprintf('  In %s>%s at %u', err.stack(i).file, err.stack(i).name, err.stack(i).line));
    end
	fprintf(2, '%s failed\n', function_name); 
	global MATLAB_RETURN_FILE; 
	MATLAB_RETURN_FILE=-1; 
end; 
exit;
