% Index the retained Cartesian pair tensors in both physical directions.
% Syntax:
%
%              [index,tensors]=cce_pair_index(spin_system)
%
% Parameters:
%
%    spin_system - spin system returned by cce_compile
%
% Outputs:
%
%    index   - sparse map from ordered physical IDs to tensor offsets
%
%    tensors - cell row; tensors{1+index(a,b)} acts on sites a,b,
%              with a shared zero tensor for absent pairs
%
% talos@spindynamics.org

function [index,tensors]=cce_pair_index(spin_system)

% Check consistency
grumble(spin_system);

% Index the reverse direction after the forward records
ids=reshape([spin_system.cce.pairs.ids],2,[]);
n=spin_system.comp.nspins;
count=size(ids,2);
index=sparse([ids(1,:) ids(2,:)],[ids(2,:) ids(1,:)],1:2*count,n,n);
forward={spin_system.cce.pairs.retained};
tensors=[{zeros(3)} forward cellfun(@transpose,forward,'UniformOutput',false)];

end

% Consistency enforcement
function grumble(spin_system)
if ~isstruct(spin_system)||~isfield(spin_system,'cce')||~isfield(spin_system.cce,'pairs')
    error('a compiled CCE spin_system is required.');
end
end


