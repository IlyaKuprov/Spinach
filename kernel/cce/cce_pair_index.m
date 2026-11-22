% Index the retained Cartesian pair tensors in both physical directions.
% Syntax:
%
%              [index,tensors]=cce_pair_index(model)
%
% Parameters:
%
%    model - CCE model returned by cce_compile
%
% Outputs:
%
%    index   - sparse map from ordered physical IDs to tensor offsets
%
%    tensors - cell row; tensors{1+index(a,b)} acts on sites a,b,
%              with a shared zero tensor for absent pairs
%
% talos@spindynamics.org

function [index,tensors]=cce_pair_index(model)

% Check consistency
grumble(model);

% Index the reverse direction after the forward records
ids=reshape([model.pairs.ids],2,[]);
n=numel(model.physical_ids);
count=size(ids,2);
index=sparse([ids(1,:) ids(2,:)],[ids(2,:) ids(1,:)],1:2*count,n,n);
forward={model.pairs.retained};
tensors=[{zeros(3)} forward cellfun(@transpose,forward,'UniformOutput',false)];

end

% Consistency enforcement
function grumble(model)
if ~isstruct(model)||~all(isfield(model,{'physical_ids','pairs'}))
    error('a compiled CCE model is required.');
end
end


