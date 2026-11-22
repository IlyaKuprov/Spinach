% Read a directed compiled Cartesian coupling for a CCE test.
% Syntax:
%
%                tensor=cce_test_read(system,ids)
%
% Parameters:
%
%    system - compiled spin system
%    ids    - ordered pair of physical spin indices
%
% Outputs:
%
%    tensor - 3-by-3 coupling in rad/s, zero for absent directions
%
% talos@spindynamics.org

function tensor=cce_test_read(system,ids)

% Check consistency
grumble(system,ids);

% Select the ordered direction without combining its reverse
list=system.inter.coupling.matrix;
row=ismember(reshape(cell2mat(list(:,1:2)),[],2),ids,'rows');
tensor=zeros(3);
if any(row)
    tensor=list{row,3};
end

end

% Consistency enforcement
function grumble(system,ids)
if ~isfield(system,'inter')||~isequal(size(ids),[1 2])
    error('a compiled spin system and an ordered ID pair are required.');
end
end


