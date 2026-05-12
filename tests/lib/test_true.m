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
% ilya.kuprov@weizmann.ac.il

function result=test_true(result,label,condition,why)

% Check the condition
if ~isscalar(condition)||~condition
    error(['FAILED: ' label ' -- ' why]);
end

% Record the pass message
result.messages{end+1}=['PASS: ' label ' -- ' why];

end
