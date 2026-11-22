% Adds coupling tensors with the same ordered endpoints. Syntax:
%
%                    pairs=sum_couplings(pairs)
%
% Parameters:
%
%    pairs - K-by-3 cell array of {from,to,tensor} rows; tensors are 3-by-3
%            Cartesian matrices in rad/s
%
% Outputs:
%
%    pairs - one sum per ordered key, sorted by to, then by from; equal
%            keys are accumulated in their original input order
%
% talos@spindynamics.org

function pairs=sum_couplings(pairs)

% Check consistency
grumble(pairs);
if isempty(pairs), return; end

% Group contributions by ordered spin pair
[keys,order]=sortrows(cell2mat(pairs(:,1:2)),[2 1]);
pairs=pairs(order,:);
starts=[1;find(any(diff(keys,1,1),2))+1];
stops=[starts(2:end)-1;size(pairs,1)];

% Sum only keys with more than one contribution
for n=find(stops>starts)'
    for k=(starts(n)+1):stops(n)
        pairs{starts(n),3}=pairs{starts(n),3}+pairs{k,3};
    end
end
pairs=pairs(starts,:);

end

% Consistency enforcement
function grumble(pairs)
if (~iscell(pairs))||(~ismatrix(pairs))||(size(pairs,2)~=3)
    error('pairs must be a three-column cell array.');
end
end


