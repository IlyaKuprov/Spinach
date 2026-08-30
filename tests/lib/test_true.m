% Adds a logical regression check with a clear message. Syntax:
%
%                    result=test_true(result,label,condition,why)
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
% A failed check is recorded in the messages and failures fields of the
% result structure rather than thrown, so that the checks that follow it
% in the calling test are still evaluated and the earlier passes are not
% lost; run_tests inspects the failures field to decide the status.
%
% ilya.kuprov@weizmann.ac.il

function result=test_true(result,label,condition,why)

% Record a failed condition and return
if ~isscalar(condition)||~condition
    result.messages{end+1}=['FAIL: ' label ' -- ' why];
    result.failures{end+1}=[label ' -- ' why];
    return
end

% Record the pass message
result.messages{end+1}=['PASS: ' label ' -- ' why];

end

