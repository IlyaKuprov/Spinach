% Adds coupling tensors with the same ordered endpoints. Syntax:
%
%                    pairs=sum_couplings(pairs)
%
% Parameters and outputs:
%
%    pairs - cell rows {from,to,tensor} in rad/s; output contains one
%            sum per key, sorted by to, then by from. Equal keys are
%            summed in their original order.

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
if (~iscell(pairs))||(size(pairs,2)~=3)
    error('pairs must be a three-column cell array.');
end
end


