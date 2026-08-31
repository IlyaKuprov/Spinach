% Runs one Spinach regression test by identifier substring. Syntax:
%
%                    result=run_test(test_id)
%
% Parameters:
%
%     test_id  - test identifier or unique substring from list_tests()
%
% Outputs:
%
%     result   - single test result structure; only returned for a
%                passing test because run_tests() throws an error
%                when the matched test fails
%
% ilya.kuprov@weizmann.ac.il

function result=run_test(test_id)

% Check consistency
grumble(test_id);

% Add the test library to the path
root_dir=fileparts(mfilename('fullpath'));
addpath(fullfile(root_dir,'lib'));

% Require a unique manifest match before running anything
manifest=test_manifest();
matched=contains({manifest.id},test_id)|...
        contains({manifest.name},test_id);
if nnz(matched)~=1
    error('test_id must match exactly one test.');
end

% Run the matching test verbosely
results=run_tests('pattern',test_id,'verbose',true,'stop_on_fail',true);
result=results;

end

% Consistency enforcement
function grumble(test_id)
if (~ischar(test_id))||isempty(test_id)||(~isrow(test_id))
    error('test_id must be a non-empty character string.');
end
end

