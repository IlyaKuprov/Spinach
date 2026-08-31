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
% run_tests and leave the remaining checks unevaluated. The record is
% retained by test_record after every check, so the catch path of
% run_tests still reports the checks completed before such an error.
%
% ilya.kuprov@weizmann.ac.il

function [result,passed]=test_true(result,label,condition,why)

% Check consistency
grumble(result,label,condition,why);

% Evaluate the condition
passed=isscalar(condition)&&(condition~=0);

% Record the outcome of the check
if passed
    result.messages{end+1}=['PASS: ' label ' -- ' why];
else
    result.messages{end+1}=['FAIL: ' label ' -- ' why];
    result.failures{end+1}=[label ' -- ' why];
end

% Retain the record for the run_tests catch path
test_record(result);

end

% Consistency enforcement
function grumble(result,label,condition,why)
if (~isstruct(result))||(~isscalar(result))||...
   (~isfield(result,'messages'))||(~isfield(result,'failures'))
    error('result must be a scalar test result structure.');
end
if (~ischar(label))||isempty(label)||(~isrow(label))
    error('label must be a non-empty character string.');
end
if (~islogical(condition))&&(~isnumeric(condition))
    error('condition must be logical or numeric.');
end
if (~ischar(why))||isempty(why)||(~isrow(why))
    error('why must be a non-empty character string.');
end
end


