% Builds a cluster Hamiltonian in rad/s. Syntax:
%
%                [H,meta]=cce_hamiltonian(spin_system,parameters,sites,spin_config)
%
% spin_system is returned by cce_compile; sites are physical bath spin IDs.
% Ordering is parameters.observed followed by sites, in the supplied order.
% spin_config is a scalar struct with optional backend (auto/direct/stock),
% exterior_ids and exterior_expectations (one real [Ix Iy Iz] row per ID).
% Exterior fields use exactly the retained pair interactions.
% meta reports IDs, dimensions, operators, backend, and roundoff residuals.
%
% Parameters:
%
%   spin_system - output of cce_compile.
%   parameters - parameters resolved by cce_options.
%   sites - physical bath IDs; [] selects the observed baseline.
%   spin_config - scalar structure with optional backend override, exterior_ids
%                 and exterior_expectations (one real Ix/Iy/Iz row per exterior ID).
%                 Exterior fields use retained pair interactions.
%
% Outputs:
%
%   H - Hermitian cluster Hamiltonian in rad/s.
%   meta - physical tensor order, dimensions, local operators, and backend
%          metadata.
%
% See also cce_compile, cce_subsystem.
%
% ilya.kuprov@weizmann.ac.il

function [H,meta]=cce_hamiltonian(spin_system,parameters,sites,spin_config)

% Check cluster membership, spin_config, and the requested backend
grumble(spin_system,parameters,sites,spin_config);

% Resolve the validated construction path
backend=parameters.cce.backend;
if isfield(spin_config,'backend')
    backend=spin_config.backend;
end
if strcmp(backend,'auto')
    if spin_system.cce.direct_supported
        backend='direct';
    else
        backend='stock';
    end
end

% Bound simultaneous dense matrices before allocating cluster operators
ids=[parameters.observed sites];
dims=spin_system.comp.mults(ids);
dimension=prod(dims);
copies=4+3*numel(ids);
if strcmp(backend,'stock')

    % Reserve every stock spherical component as a dense worst case
    ranks=max([2 cellfun(@numel,spin_system.inter.giant.coeff).']);
    copies=copies+sum((2*(1:ranks)+1).^2);

end
if ~isfinite(dimension)||16*dimension^2*copies>spin_system.tols.cce_mem
    error('cluster Hamiltonian exceeds sys.tols.cce_mem.');
end

% Resolve both directions of each retained pair
[pair_index,tensors]=cce_pair_index(spin_system);

% Construct angular momentum operators in canonical product ordering
ops=cell(numel(ids),3);
for n=1:numel(ids)

    % Build local raising, lowering, and longitudinal matrices
    spin=(dims(n)-1)/2;
    magnetic=spin:-1:-spin;
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
    H=sparse(dimension,dimension);
    sum_norm=0;
    operations=0;
    for n=1:numel(ids)

        % Select Zeeman components and rotating-frame carrier convention
        id=ids(n);
        tensor=spin_system.cce.zeeman{id};
        strength=spin_system.inter.zeeman.strength{id};
        field=zeros(1,3);
        isotropic=trace(spin_system.inter.zeeman.matrix{id})/3;
        anisotropic=tensor-isotropic*eye(3);
        [~,rank_one,rank_two]=mat2sphten(spin_system.inter.zeeman.matrix{id});
        if norm(rank_one,2)<=spin_system.tols.liouv_zero&&...
           norm(rank_two,2)<=spin_system.tols.liouv_zero
            anisotropic=zeros(3);
        end
        switch strength
            case 'full'
                if abs(isotropic)<=spin_system.tols.liouv_zero
                    isotropic=0;
                end
                field=anisotropic(3,:)+[0 0 isotropic];
            case 'secular'
                isotropic=isotropic-spin_system.inter.basefrqs(id);
                if abs(isotropic)<=spin_system.tols.liouv_zero
                    isotropic=0;
                end
                field(3)=anisotropic(3,3)+isotropic;
            case 'ignore'

                % No local Zeeman term retained
            otherwise
                error('unvalidated direct Zeeman retention.');
        end
        for axis=1:3
            term=field(axis)*ops{n,axis};
            H=H+term;
            sum_norm=sum_norm+norm(term,'fro');
            operations=operations+1;
        end

        % Apply the stock on-site quadratic retention convention
        tensor=spin_system.cce.quadratic{id};
        strength=spin_system.cce.quadratic_strength{id};
        switch strength
            case 'strong'

                % Keep the full symmetric traceless tensor
            case 'secular'
                tensor=diag([-tensor(3,3)/2 -tensor(3,3)/2 tensor(3,3)]);
            case 'zz'
                if norm(tensor,2)>2*pi*spin_system.tols.inter_cutoff
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
                term=tensor(axis,other)*(ops{n,axis}*ops{n,other});
                H=H+term;
                sum_norm=sum_norm+norm(term,'fro');
                operations=operations+dimension;
            end
        end

        % Sum each physical pair once, using combined ordered-cell tensors
        for k=n+1:numel(ids)
            tensor=tensors{1+pair_index(id,ids(k))};
            for axis=1:3
                for other=1:3
                    term=tensor(axis,other)*(ops{n,axis}*ops{k,other});
                    H=H+term;
                    sum_norm=sum_norm+norm(term,'fro');
                    operations=operations+dimension;
                end
            end
        end

    end
else

    % Remove an inherited basis before slicing immutable physical metadata
    system=spin_system;
    if isfield(system,'bas')
        system=rmfield(system,'bas');
    end
    system=kill_spin(system,setdiff(1:spin_system.comp.nspins,ids));
    basis_parameters.formalism='zeeman-hilb';
    basis_parameters.approximation={'none'};
    system=basis(system,basis_parameters);
    system=assume(system,spin_system.cce.assumptions{:});
    [H,Q]=hamiltonian(system);
    H=H+orientation(Q,parameters.orientation);

    % Reorder the stock ascending physical indices into observed-first order
    sorted_ids=sort(ids);
    sorted_dims=spin_system.comp.mults(sorted_ids);
    [~,order]=ismember(ids,sorted_ids);
    if numel(ids)>1
        indices=reshape(1:dimension,fliplr(sorted_dims));
        indices=permute(indices,numel(ids)+1-fliplr(order));
        indices=indices(:);
        H=H(indices,indices);
    end
    sum_norm=norm(H,'fro');
    operations=dimension*numel(ids)^2;

end

% Add static exterior mean fields after the common ordering is established
if isfield(spin_config,'exterior_ids')
    for n=1:numel(ids)
        for k=1:numel(spin_config.exterior_ids)
            field=tensors{1+pair_index(ids(n),spin_config.exterior_ids(k))}*...
                  spin_config.exterior_expectations(k,:).';
            for axis=1:3
                term=field(axis)*ops{n,axis};
                H=H+term;
                sum_norm=sum_norm+norm(term,'fro');
                operations=operations+dimension;
            end
        end
    end
end

% Report the unsymmetrised residual and a dimension/operation roundoff bound
meta.antihermitian_residual=norm(H-H','fro');
roundoff=operations*eps/(1-operations*eps);
meta.roundoff_bound=roundoff*max(sum_norm,realmin)+...
                    dimension*spin_system.cce.canonical_roundoff;
if meta.antihermitian_residual>2*meta.roundoff_bound
    error('Hamiltonian has a non-roundoff anti-Hermitian residual.');
end
H=(H+H')/2;

% Return ordering, operators, backend, and approximation provenance
meta.physical_ids=ids;
meta.dimensions=dims;
meta.dimension=dimension;
meta.observed_dimension=prod(spin_system.comp.mults(parameters.observed));
meta.bath_dimension=prod(spin_system.comp.mults(sites));
meta.operators=ops;
meta.backend=backend;
meta.fallback_reason=spin_system.cce.fallback_reason;
meta.source_hash=spin_system.cce.source_hash;
meta.omitted_pairs=spin_system.cce.omitted_pairs;
meta.observable_budget=spin_system.tols.cce_atol+spin_system.tols.cce_rtol*norm(H,'fro');

end

% Validate structure, IDs, and explicitly supplied exterior expectations
function grumble(spin_system,parameters,sites,spin_config)
if ~isstruct(spin_system)||~isfield(spin_system,'cce')||...
   ~all(isfield(spin_system.cce,{'source_hash','pairs'}))
    error('model must be returned by cce_compile.');
end
if ~isnumeric(sites)||~isreal(sites)||(~isempty(sites)&&~isrow(sites))||...
   any(~isfinite(sites))||any(sites~=fix(sites))||...
   numel(unique(sites))~=numel(sites)||~all(ismember(sites,parameters.cce.bath_spins))
    error('sites must be a row of unique spin_system bath spin IDs.');
end
if ~isstruct(spin_config)||~isscalar(spin_config)||...
   ~all(ismember(fieldnames(spin_config),{'backend','exterior_ids','exterior_expectations'}))
    error('spin_config contains unknown fields or is not a scalar structure.');
end
if isfield(spin_config,'backend')&&~ismember(spin_config.backend,{'auto','direct','stock'})
    error('spin_config.backend must be auto, direct, or stock.');
end
if ~spin_system.cce.direct_supported&&...
   ((isfield(spin_config,'backend')&&strcmp(spin_config.backend,'direct'))||...
    (~isfield(spin_config,'backend')&&strcmp(parameters.cce.backend,'direct')))
    error('direct Hamiltonian is not validated for %s.',spin_system.cce.fallback_reason);
end
if xor(isfield(spin_config,'exterior_ids'),isfield(spin_config,'exterior_expectations'))
    error('exterior_ids and exterior_expectations must be supplied together.');
end
if isfield(spin_config,'exterior_ids')
    ids=spin_config.exterior_ids;
    values=spin_config.exterior_expectations;
    if ~isnumeric(ids)||~isreal(ids)||(~isempty(ids)&&~isrow(ids))||...
       any(~isfinite(ids))||numel(unique(ids))~=numel(ids)||...
       ~all(ismember(ids,setdiff(1:spin_system.comp.nspins,[parameters.observed sites])))||...
       ~isnumeric(values)||~isreal(values)||~isequal(size(values),[numel(ids) 3])||...
       any(~isfinite(values),'all')
        error('invalid exterior IDs or expectation rows.');
    end
    if any(vecnorm(values,2,2)>(spin_system.comp.mults(ids).'-1)/2+eps)
        error('exterior expectation exceeds the physical spin length.');
    end
    if ~spin_system.cce.direct_supported&&~isempty(ids)
        error('exterior mean fields require validated pair-retention semantics.');
    end
end
end


