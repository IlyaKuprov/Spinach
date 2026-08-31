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
% A test is reported as a failure when it records at least one failed
% check in the failures field of its result structure, in which case the
% error field carries the recorded failures; the messages accumulated
% before the failure are preserved. A test that throws a MATLAB error is
% also reported as a failure, and the suite continues with the next test;
% the checks it had recorded before the error are recovered from
% test_record and reported alongside the error message.
%
% ilya.kuprov@weizmann.ac.il

function results=run_tests(varargin)

% Check consistency
grumble(varargin);

% Add the test library to the path
root_dir=fileparts(mfilename('fullpath'));
addpath(root_dir);
addpath(fullfile(root_dir,'lib'));
addpath(genpath(fullfile(root_dir,'kernel')));
addpath(genpath(fullfile(root_dir,'interfaces')));

% Add the Spinach production directories to the path
spinach_root=fileparts(root_dir);
addpath(genpath(fullfile(spinach_root,'etc')));
addpath(genpath(fullfile(spinach_root,'experiments')));
addpath(genpath(fullfile(spinach_root,'interfaces')));
addpath(genpath(fullfile(spinach_root,'kernel')));

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
               'messages',{},'failures',{},'error',{});

% Run the tests
for n=1:numel(manifest)

    % Reset the clock and the record accumulator
    tic;
    test_record(new_test_result(manifest(n).id,manifest(n).name,''));

    try
        result=feval(manifest(n).function);
        result.elapsed=toc;
        result.error=strjoin(result.failures,'; ');
        if isempty(result.failures)
            result.status='PASS';
        else
            result.status='FAIL';
        end
    catch err
        result=test_record([]);
        result.status='FAIL';
        result.elapsed=toc;
        result.failures{end+1}=err.message;
        result.error=strjoin(result.failures,'; ');
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

% Consistency enforcement
function grumble(option_list)
if mod(numel(option_list),2)~=0
    error('options must be supplied as name-value pairs.');
end
for n=1:2:numel(option_list)
    if (~ischar(option_list{n}))||isempty(option_list{n})||(~isrow(option_list{n}))
        error('option names must be non-empty character strings.');
    end
end
end


