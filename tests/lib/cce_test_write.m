% Install a directed Cartesian tensor in a compiled CCE test fixture.
% Syntax:
%
%                system=cce_test_write(system,ids,tensor)
%
% Parameters:
%
%    system - compiled spin system before assumption assignment
%    ids    - ordered pair of physical spin indices
%    tensor - 3-by-3 coupling in rad/s; zero deletes the direction
%
% Outputs:
%
%    system - modified fixture with canonical ordered coupling records
%
% talos@spindynamics.org

function system=cce_test_write(system,ids,tensor)

% Check consistency
grumble(system,ids,tensor);

% Replace the old record and restore compiled destination-major ordering
list=system.inter.coupling.matrix;
keys=reshape(cell2mat(list(:,1:2)),[],2);
list(ismember(keys,ids,'rows'),:)=[];
if any(tensor,'all')
    list(end+1,:)={ids(1),ids(2),tensor};
end
system.inter.coupling.matrix=sortrows(list,[2 1]);

end

% Consistency enforcement
function grumble(system,ids,tensor)
if ~isfield(system,'inter')||~isequal(size(ids),[1 2])||~isequal(size(tensor),[3 3])
    error('a compiled system, ordered ID pair, and Cartesian tensor are required.');
end
end


