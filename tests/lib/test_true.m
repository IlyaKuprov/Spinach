% Adds a logical regression check with a clear message. Syntax:
%
%              [result,passed]=test_true(result,label,condition,why)
%
% Parameters:
%
%     result     - test result structure
%
%     label      - check label
%
%     condition  - logical pass/fail condition
%
%     why        - explanation of the right answer
%
% Outputs:
%
%     result     - updated test result structure
%
%     passed     - true when the condition held
%
% A failed check is recorded in the messages and failures fields of the
% result structure rather than thrown, so that the checks that follow it
% in the calling test are still evaluated and the earlier passes are not
% lost; run_tests inspects the failures field to decide the status. Where
% the check is a precondition for the statements that follow it, the
% caller must guard those statements with the passed flag, because a
% secondary Matlab error would send the test into the catch path of
% run_tests and discard the informative message recorded here.
%
% ilya.kuprov@weizmann.ac.il

function [result,passed]=test_true(result,label,condition,why)

% Evaluate the condition
passed=isscalar(condition)&&(condition~=0);

% Record a failed condition and return
if ~passed
    result.messages{end+1}=['FAIL: ' label ' -- ' why];
    result.failures{end+1}=[label ' -- ' why];
    return
end

% Record the pass message
result.messages{end+1}=['PASS: ' label ' -- ' why];

end

