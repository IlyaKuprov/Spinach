% Keeps the test result under construction available to the harness. Syntax:
%
%                        record=test_record(record)
%
% Parameters:
%
%     record  - test result structure to retain, or an empty array to
%               retrieve the structure retained by the previous call
%
% Outputs:
%
%     record  - the structure just retained, or the structure retained
%               by the previous call when the input is empty
%
% Test functions pass their result structure by value, so a Matlab error
% inside a test destroys every check it had already recorded. run_tests
% retains a blank record before each test, test_true and test_close
% retain the updated record after each check, and the catch block of
% run_tests retrieves what survived and reports it alongside the error.
%
% ilya.kuprov@weizmann.ac.il

function record=test_record(record)

% Retain the record between calls
persistent latest

% Check consistency
grumble(record);

% Retrieve the retained record or retain the supplied one
if isempty(record)
    record=latest;
else
    latest=record;
end

end

% Consistency enforcement
function grumble(record)
if (~isempty(record))&&((~isstruct(record))||(~isscalar(record)))
    error('record must be a scalar structure or an empty array.');
end
end


