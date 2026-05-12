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
%     result   - test result structure
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
result.error='';

end
