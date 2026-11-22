% Restricts compiled pair records to retained particles. Syntax:
%
%                 pairs=pair_subset(pairs,spin_map)
%
% Parameters:
%
%    pairs    - K-by-3 cell array of ordered {from,to,value} records
%
%    spin_map - column vector mapping old particle indices to new
%               indices; zero means removed, and positive entries
%               are consecutive in ascending order
%
% Outputs:
%
%    pairs    - retained records with reindexed endpoints, in the
%               original record order
%
% talos@spindynamics.org

function pairs=pair_subset(pairs,spin_map)

% Check consistency
grumble(pairs,spin_map);
if isempty(pairs), return; end

% Drop records touching removed particles and remap the remaining keys
keys=cell2mat(pairs(:,1:2));
keys=reshape(spin_map(keys),size(keys));
keep=all(keys>0,2);
pairs=pairs(keep,:);
pairs(:,1:2)=num2cell(keys(keep,:));

end

% Consistency enforcement
function grumble(pairs,spin_map)
if (~iscell(pairs))||(~ismatrix(pairs))||(size(pairs,2)~=3)
    error('pairs must be a three-column cell array.');
end
if (~isnumeric(spin_map))||(~isreal(spin_map))||(~iscolumn(spin_map))
    error('spin_map must be a real numeric column vector.');
end
end


