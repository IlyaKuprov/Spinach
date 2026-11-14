% Builds a cluster Hamiltonian in rad/s. Syntax:
%
%                [H,meta]=cce_hamiltonian(model,sites,context)
%
% model is returned by cce_prepare; sites are physical bath spin IDs.
% Ordering is model.central_ids followed by sites, in the supplied order.
% context is a scalar struct with optional backend (auto/direct/stock),
% exterior_ids and exterior_expectations (one real [Ix Iy Iz] row per ID).
% Exterior fields use exactly the model's retained pair interactions.
% meta reports IDs, dimensions, operators, backend and roundoff residuals.
%
% Inputs:
%   model - output of cce_prepare.
%   sites - physical bath IDs; [] selects the central baseline.
%   context - scalar structure with optional backend override, exterior_ids
%   and exterior_expectations (one real Ix/Iy/Iz row per exterior ID).
%   Exterior fields use retained model pair interactions.
%
% Outputs:
%   H - Hermitian cluster Hamiltonian in rad/s.
%   meta - physical tensor order, dimensions, local operators and backend
%   metadata.
%
% See also cce_prepare, cce_states.
%
% ilya.kuprov@weizmann.ac.il

function [H,meta]=cce_hamiltonian(model,sites,context)

% Check cluster membership, context and the requested backend
grumble(model,sites,context);

% Resolve the validated construction path
backend=model.backend;
if isfield(context,'backend'), backend=context.backend; end
if strcmp(backend,'auto')
    if model.direct_supported, backend='direct'; else, backend='stock'; end
end
if strcmp(backend,'direct')&&~model.direct_supported
    error('direct Hamiltonian is not validated for %s.',model.fallback_reason);
end

% Bound simultaneous dense matrices before allocating cluster operators
ids=[model.central_ids sites]; dims=model.multiplicities(ids);
dimension=prod(dims); copies=4+3*numel(ids);
if strcmp(backend,'stock')

    % Reserve every stock spherical component as a dense worst case
    ranks=max([2 cellfun(@numel,model.spin_system.inter.giant.coeff).']);
    copies=copies+sum((2*(1:ranks)+1).^2);
end
if ~isfinite(dimension)||16*dimension^2*copies>model.numerics.max_memory_bytes
    error('cluster Hamiltonian exceeds numerics.max_memory_bytes.');
end

% Construct angular momentum operators in canonical product ordering
ops=cell(numel(ids),3);
for n=1:numel(ids)

    % Build local raising, lowering, and longitudinal matrices
    spin=(dims(n)-1)/2; magnetic=spin:-1:-spin;
    raising=diag(sqrt(spin*(spin+1)-magnetic(2:end).*(magnetic(2:end)+1)),1);
    local={(raising+raising')/2,(raising-raising')/(2i),diag(magnetic)};

    % Embed one local matrix at a time into the cluster product space
    for axis=1:3
        ops{n,axis}=kron(kron(speye(prod(dims(1:n-1))),...
                            sparse(local{axis})),speye(prod(dims(n+1:end))));
    end
end

% Build either the validated Cartesian Hamiltonian or the stock reference
if strcmp(backend,'direct')

    % Accumulate all retained local and unordered pair terms
    H=sparse(dimension,dimension); sum_norm=0; operations=0;
    for n=1:numel(ids)

        % Select Zeeman components and rotating-frame carrier convention
        id=ids(n); tensor=model.zeeman{id};
        strength=model.spin_system.inter.zeeman.strength{id};
        field=zeros(1,3); isotropic=trace(model.spin_system.inter.zeeman.matrix{id})/3;
        anisotropic=tensor-isotropic*eye(3);
        [~,rank_one,rank_two]=mat2sphten(model.spin_system.inter.zeeman.matrix{id});
        if norm(rank_one,2)<=model.spin_system.tols.liouv_zero&&...
           norm(rank_two,2)<=model.spin_system.tols.liouv_zero
            anisotropic=zeros(3);
        end
        switch strength
            case 'full'
                if abs(isotropic)<=model.spin_system.tols.liouv_zero, isotropic=0; end
                field=anisotropic(3,:)+[0 0 isotropic];
            case 'secular'
                isotropic=isotropic-model.spin_system.inter.basefrqs(id);
                if abs(isotropic)<=model.spin_system.tols.liouv_zero, isotropic=0; end
                field(3)=anisotropic(3,3)+isotropic;
            case 'ignore'

                % No local Zeeman term retained
            otherwise
                error('unvalidated direct Zeeman retention.');
        end
        for axis=1:3
            term=field(axis)*ops{n,axis}; H=H+term;
            sum_norm=sum_norm+norm(term,'fro'); operations=operations+1;
        end

        % Apply the stock on-site quadratic retention convention
        tensor=model.quadratic{id};
        strength=model.spin_system.inter.coupling.strength{id,id};
        switch strength
            case 'strong'

                % Keep the full symmetric traceless tensor
            case 'secular'
                tensor=diag([-tensor(3,3)/2 -tensor(3,3)/2 tensor(3,3)]);
            case 'zz'
                if norm(tensor,2)>2*pi*model.spin_system.tols.inter_cutoff
                    error('zz on-site quadratic retention is undefined by Spinach.');
                end
                tensor=zeros(3);
            case 'ignore'
                tensor=zeros(3);
            otherwise
                error('unvalidated direct quadratic retention.');
        end
        for axis=1:3
            for other=1:3
                term=tensor(axis,other)*(ops{n,axis}*ops{n,other}); H=H+term;
                sum_norm=sum_norm+norm(term,'fro'); operations=operations+dimension;
            end
        end

        % Sum each physical pair once, using combined ordered-cell tensors
        for k=n+1:numel(ids)
            tensor=model.retained{id,ids(k)};
            for axis=1:3
                for other=1:3
                    term=tensor(axis,other)*(ops{n,axis}*ops{k,other}); H=H+term;
                    sum_norm=sum_norm+norm(term,'fro'); operations=operations+dimension;
                end
            end
        end
    end
else

    % Remove an inherited basis before slicing immutable physical metadata
    system=model.spin_system;
    if isfield(system,'bas'), system=rmfield(system,'bas'); end
    system=kill_spin(system,setdiff(model.physical_ids,ids));
    basis_parameters.formalism='zeeman-hilb'; basis_parameters.approximation={'none'};
    system=basis(system,basis_parameters);
    system=assume(system,model.assumptions{:});
    [H,Q]=hamiltonian(system); H=H+orientation(Q,model.orientation);

    % Reorder the stock ascending physical indices into central-first order
    sorted_ids=sort(ids); sorted_dims=model.multiplicities(sorted_ids);
    [~,order]=ismember(ids,sorted_ids);
    if numel(ids)>1
        indices=reshape(1:dimension,fliplr(sorted_dims));
        indices=permute(indices,numel(ids)+1-fliplr(order));
        indices=indices(:); H=H(indices,indices);
    end
    sum_norm=norm(H,'fro'); operations=dimension*numel(ids)^2;
end

% Add static exterior mean fields after the common ordering is established
if isfield(context,'exterior_ids')
    if ~model.direct_supported&&~isempty(context.exterior_ids)
        error('exterior mean fields require validated pair-retention semantics.');
    end
    for n=1:numel(ids)
        for k=1:numel(context.exterior_ids)
            field=model.retained{ids(n),context.exterior_ids(k)}*...
                  context.exterior_expectations(k,:).';
            for axis=1:3
                term=field(axis)*ops{n,axis}; H=H+term;
                sum_norm=sum_norm+norm(term,'fro'); operations=operations+dimension;
            end
        end
    end
end

% Report the unsymmetrized residual and a dimension/operation roundoff bound
meta.antihermitian_residual=norm(H-H','fro');
roundoff=operations*eps/(1-operations*eps);
meta.roundoff_bound=roundoff*max(sum_norm,realmin)+...
                    dimension*model.canonical_roundoff;
if meta.antihermitian_residual>2*meta.roundoff_bound
    error('Hamiltonian has a non-roundoff anti-Hermitian residual.');
end
H=(H+H')/2;

% Return ordering, operators, backend and approximation provenance
meta.physical_ids=ids; meta.dimensions=dims; meta.dimension=dimension;
meta.central_dimension=prod(model.multiplicities(model.central_ids));
meta.bath_dimension=prod(model.multiplicities(sites)); meta.operators=ops;
meta.backend=backend; meta.fallback_reason=model.fallback_reason;
meta.source_hash=model.source_hash; meta.omitted_pairs=model.omitted_pairs;
meta.observable_budget=model.numerics.atol+model.numerics.rtol*norm(H,'fro');

end

% Validate structure, IDs, and explicitly supplied exterior expectations
function grumble(model,sites,context)
if ~isstruct(model)||~all(isfield(model,{'source_hash','central_ids',...
   'bath_ids','multiplicities','retained','spin_system','numerics','backend'}))
    error('model must be returned by cce_prepare.');
end
if ~isnumeric(sites)||~isreal(sites)||(~isempty(sites)&&~isrow(sites))||...
   any(~isfinite(sites))||any(sites~=fix(sites))||...
   numel(unique(sites))~=numel(sites)||~all(ismember(sites,model.bath_ids))
    error('sites must be a row of unique model bath spin IDs.');
end
if ~isstruct(context)||~isscalar(context)||...
   ~all(ismember(fieldnames(context),{'backend','exterior_ids','exterior_expectations'}))
    error('context contains unknown fields or is not a scalar structure.');
end
if isfield(context,'backend')&&~ismember(context.backend,{'auto','direct','stock'})
    error('context.backend must be auto, direct, or stock.');
end
if xor(isfield(context,'exterior_ids'),isfield(context,'exterior_expectations'))
    error('exterior_ids and exterior_expectations must be supplied together.');
end
if isfield(context,'exterior_ids')
    ids=context.exterior_ids; values=context.exterior_expectations;
    if ~isnumeric(ids)||~isreal(ids)||(~isempty(ids)&&~isrow(ids))||...
       any(~isfinite(ids))||numel(unique(ids))~=numel(ids)||...
       ~all(ismember(ids,setdiff(model.physical_ids,[model.central_ids sites])))||...
       ~isnumeric(values)||~isreal(values)||~isequal(size(values),[numel(ids) 3])||...
       any(~isfinite(values),'all')
        error('invalid exterior IDs or expectation rows.');
    end
    if any(vecnorm(values,2,2)>(model.multiplicities(ids).'-1)/2+eps)
        error('exterior expectation exceeds the physical spin length.');
    end
end
end


