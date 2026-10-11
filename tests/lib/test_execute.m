% Executes one regression test with captured output and isolated graphics.
% Syntax:
%
%                  [result,transcript]=test_execute(entry)
%
% Parameters:
%
%     entry      - scalar manifest entry with id, name, and function
%
% Outputs:
%
%     result     - standard test result, including recorded and thrown
%                  failures and elapsed time
%
%     transcript - captured console output, including output before errors

function [result,transcript]=test_execute(entry)

% Check consistency
grumble(entry);

% Isolate graphics and retain partial results when a test throws
cleaner=test_graphics(); %#ok<NASGU>
clock=tic;
test_record(new_test_result(entry.id,entry.name,''));

% Catch inside evalc so output preceding an exception is retained
transcript=evalc(['try; result=feval(entry.function); '...
                 'catch err; result=test_record([]); '...
                 'result.failures{end+1}=err.message; end']);

% Clean graphics before returning the outcome to the suite
clear('cleaner');
result.elapsed=toc(clock);
result.error=strjoin(result.failures,'; ');
if isempty(result.failures)
    result.status='PASS';
else
    result.status='FAIL';
end

end

% Consistency enforcement
function grumble(entry)
if (~isstruct(entry))||(~isscalar(entry))||...
   (~all(isfield(entry,{'id','name','function'})))
    error('entry must be a scalar test manifest structure.');
end
if (~ischar(entry.id))||isempty(entry.id)||(~isrow(entry.id))||...
   (~ischar(entry.name))||isempty(entry.name)||(~isrow(entry.name))
    error('test id and name must be non-empty character row vectors.');
end
if (~isa(entry.function,'function_handle'))&&...
   ((~ischar(entry.function))||isempty(entry.function)||(~isrow(entry.function)))
    error('test function must be a function handle or a character row vector.');
end
end


