% Creates a regression test result structure. Syntax:
%
%                    result=new_test_result(id,name,purpose)
%
% Parameters:
%
%     id       - stable test identifier
%
%     name     - short human-readable test name
%
%     purpose  - one-sentence purpose statement
%
% Outputs:
%
%     result   - test result structure, in which the messages field
%                accumulates one line per check and the failures field
%                accumulates the details of the checks that did not
%                pass; an empty failures field means the test passed
%
% ilya.kuprov@weizmann.ac.il

function result=new_test_result(id,name,purpose)

% Check consistency
grumble(id,name,purpose);

% Build the result structure
result.id=id;
result.name=name;
result.purpose=purpose;
result.status='RUNNING';
result.elapsed=0;
result.messages={};
result.failures={};
result.error='';

end

% Consistency enforcement
function grumble(id,name,purpose)
if (~ischar(id))||isempty(id)||(~isrow(id))
    error('id must be a non-empty character string.');
end
if (~ischar(name))||isempty(name)||(~isrow(name))
    error('name must be a non-empty character string.');
end
if (~ischar(purpose))||((~isempty(purpose))&&(~isrow(purpose)))
    error('purpose must be a character string.');
end
end


