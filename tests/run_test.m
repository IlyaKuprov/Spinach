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
%     result   - single test result structure
%
% ilya.kuprov@weizmann.ac.il

function result=run_test(test_id)

% Run one matching test verbosely
results=run_tests('pattern',test_id,'verbose',true,'stop_on_fail',true);
if numel(results)~=1
    error('test_id must match exactly one test.');
end
result=results;

end
