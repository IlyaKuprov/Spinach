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

