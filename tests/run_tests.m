% Runs the Spinach regression test suite. Syntax:
%
%                    results=run_tests(varargin)
%
% Parameters:
%
%     varargin  - name-value options: 'pattern', 'verbose', and
%                 'stop_on_fail'
%
% Outputs:
%
%     results   - structure array with test outcomes and messages
%
% ilya.kuprov@weizmann.ac.il

function results=run_tests(varargin)

% Add the test library to the path
root_dir=fileparts(mfilename('fullpath'));
addpath(root_dir);
addpath(fullfile(root_dir,'lib'));
addpath(genpath(fullfile(root_dir,'kernel')));

% Parse options
options=test_options(varargin{:});

% Get the manifest
manifest=test_manifest();

% Apply substring filter
if ~isempty(options.pattern)
    keep=contains({manifest.id},options.pattern)|...
         contains({manifest.name},options.pattern);
    manifest=manifest(keep);
end

% Preallocate result array
results=struct('id',{},'name',{},'purpose',{},'status',{},'elapsed',{},...
               'messages',{},'error',{});

% Run the tests
for n=1:numel(manifest)
    tic;
    try
        result=feval(manifest(n).function);
        result.status='PASS';
        result.elapsed=toc;
        result.error='';
    catch err
        result.id=manifest(n).id;
        result.name=manifest(n).name;
        result.purpose='';
        result.status='FAIL';
        result.elapsed=toc;
        result.messages={};
        result.error=err.message;
    end
    results(end+1)=result; %#ok<AGROW>
    if options.verbose
        fprintf('%s\t%s\t%.3f s\n',result.status,result.id,result.elapsed);
        for k=1:numel(result.messages)
            fprintf('    %s\n',result.messages{k});
        end
        if strcmp(result.status,'FAIL')
            fprintf('    ERROR: %s\n',result.error);
        end
    end
    if options.stop_on_fail&&strcmp(result.status,'FAIL')
        break;
    end
end

% Summarise outcomes
n_pass=nnz(strcmp({results.status},'PASS'));
n_fail=nnz(strcmp({results.status},'FAIL'));
fprintf('Spinach regression tests: %d passed, %d failed.\n',n_pass,n_fail);
if n_fail>0
    failed=results(strcmp({results.status},'FAIL'));
    for n=1:numel(failed)
        fprintf('FAIL\t%s\t%s\n',failed(n).id,failed(n).error);
    end
    error('Spinach regression tests failed.');
end

end
