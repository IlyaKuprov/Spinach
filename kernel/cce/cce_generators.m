% Construct retained-cluster generators in Spinach angular-frequency units.
% Syntax:
%
%   [H,R,K,parameters,cache]=cce_generators(cluster_system,parameters,cache)
%
% cluster_system is cce_subsystem output; parameters is a resolved native
% experiment structure. Conventional generators obey L=H+1i*R+1i*K and
% exp(-1i*L*t). R is the GKSL dissipator in s^-1; K is zero. Dissipation
% requires zeeman-liouv, except explicitly tagged projected branches. A
% projected dissipator acts on an alpha/beta coherence block, not a density.
% Its swapped branch is swap_matrix*conj(dissipator)*swap_matrix'.
% parameters.bath_jumps uses support, operator, and rate; supported spin
% descriptors are L+, L-, Lx, Ly, and Lz, or explicit ordered matrices.
% Zero-rate jumps are omitted before closure checks; active split supports
% and exterior jumps driving a frozen mean field are refused.
% Projected H has kind='cce_projected', left/right Hamiltonians and a
% dissipator acting on column-vectorised bath coherence. R and K are []:
% the tagged object, not numeric addition, goes to cce_evolution. Swap
% left/right and reverse the projected dissipator with swap_matrix for
% a refocusing pulse. These tags are never stock pulse-sequence inputs.
% parameters.cce.backend selects auto/direct/stock. The returned
% parameters.cce_aux contains ordering, Hamiltonian, and rounding metadata.
% Master-equation cce_aux.generator retains the physical Hamiltonians,
% spin operators, and included jump pairs for completed raw-task evidence.
% cache is [] or the separate exact generator cache returned here.
%
% Parameters:
%
%   cluster_system - retained subsystem returned by cce_subsystem
%   parameters     - resolved native experiment and numerical options
%   cache          - [] or the exact generator cache returned previously
%
% Outputs:
%
%   H,R,K      - retained generators or tagged projected branches
%   parameters - input options with generator diagnostics in cce_aux
%   cache      - exact-input generator cache with hit and miss counts
%
% ilya.kuprov@weizmann.ac.il

function [H,R,K,parameters,cache]=cce_generators(cluster_system,parameters,cache)

% Validate ownership before allocating retained-space matrices
grumble(cluster_system,parameters,cache);
model=cluster_system.model;
ids=cluster_system.physical_ids;
d=cluster_system.dimension;
db=cluster_system.bath_dimension;
projected=strcmp(parameters.cce.approximation,'projected');
jumps=[];
if isfield(parameters,'bath_jumps')
    jumps=parameters.bath_jumps;
end

% Exact input comparison prevents stale orientation or configuration reuse
if isfield(parameters,'cce_aux')
    parameters=rmfield(parameters,'cce_aux');
end
signature=struct('cluster_system',cluster_system,'parameters',parameters);
if ~isempty(cache)&&isequaln(cache.signature,signature)
    H=cache.H;
    R=cache.R;
    K=cache.K;
    parameters=cache.parameters;
    cache.hits=cache.hits+1;
    return
end

% Count resident jump operators before Hamiltonian or dissipator allocation
if ~isempty(jumps)||(projected&&strcmp(cluster_system.formalism,'zeeman-liouv'))
    if projected
        n=db;
    else
        n=d;
    end
    included=0;
    for k=1:numel(jumps)
        if all(ismember(jumps(k).support,ids))
            included=included+1;
        end
    end
    memory_guard(parameters,16*((12+3*numel(ids))*d^2+8*n^4+2*included*d^2));
end

% Preserve the canonical direct and independent stock construction paths
[ham,meta]=matrix(cluster_system,parameters);
parameters.cce_aux=meta;
parameters.cce_aux.generator_convention='H+iR+iK';
parameters.cce_aux.hamiltonian_units='rad/s';
parameters.cce_aux.jump_rate_units='s^-1';
if projected

    % Conditional branches use the declared observed transition, never a guess
    if ~isfield(model,'transition')
        error('cce_generators:transition','Projected evolution needs a transition.');
    end
    transition=model.transition;
    va=kron(transition.alpha,eye(db));
    vb=kron(transition.beta,eye(db));
    left=va'*ham*va;
    right=vb'*ham*vb;
    leakage=max(norm(ham*va-va*left,'fro'),norm(ham*vb-vb*right,'fro'));
    parameters.cce_aux.projection_leakage=leakage;
    if ~isempty(jumps)&&leakage>d*eps*norm(ham,'fro')
        error('cce:population_transfer','Projected master dynamics requires invariant observed manifolds.');
    end
    second=parameters.cce.second_order.enabled;
    if second
        if ~isempty(jumps)
            error('cce:capability','Dissipative virtual-state corrections are unsupported.');
        end

        % Resolve isolated observed eigenvectors against the compiled transition
        configuration=rmfield(cluster_system.configuration,{'exterior_ids','exterior_expectations'});
        [observed_system,~]=cce_subsystem(model,[],configuration,[]);
        observed_system.configuration.exterior_ids=[];
        observed_system.configuration.exterior_expectations=zeros(0,3);
        [observed_hamiltonian,~]=matrix(observed_system,parameters);
        [vectors,energies]=eig(full((observed_hamiltonian+observed_hamiltonian')/2),'vector');
        [energies,order]=sort(real(energies));
        vectors=vectors(:,order);
        [~,selected(1)]=max(abs(vectors'*model.transition.alpha).^2);
        [~,selected(2)]=max(abs(vectors'*model.transition.beta).^2);
        second=parameters.cce.second_order;
        observed=model.observed_ids;
        ops=cell(1,3*numel(observed));
        fields=cell(size(ops));
        axes={'x','y','z'};
        for c=1:numel(observed)
            for a=1:3
                slot=3*(c-1)+a;
                ops{slot}=embed(pauli_axis(model.multiplicities(observed(c)),axes{a}),...
                               c,model.multiplicities(observed));
                field=zeros(cluster_system.bath_dimension);
                for b=1:numel(cluster_system.sites)
                    tensor=model.retained{observed(c),cluster_system.sites(b)};
                    for axis=1:3
                        field=field+tensor(a,axis)*embed(...
                            pauli_axis(cluster_system.multiplicities(numel(observed)+b),axes{axis}),...
                            b,cluster_system.multiplicities(numel(observed)+1:end));
                    end
                end
                for b=1:numel(cluster_system.configuration.exterior_ids)
                    field=field+eye(cluster_system.bath_dimension)*(...
                        model.retained{observed(c),cluster_system.configuration.exterior_ids(b)}(a,:)*...
                        cluster_system.configuration.exterior_expectations(b,:).');
                end
                fields{slot}=field;
            end
        end

        % Add the declared virtual-state correction to each conditional branch
        for branch=1:2
            state=vectors(:,selected(branch));
            gaps=energies(selected(branch))-energies;
            outside=setdiff(1:numel(energies),selected(branch));
            if any(abs(gaps(outside))<=second.gap_tolerance)
                error('Virtual-state energy gap is below the explicit tolerance.');
            end
            names={'alpha','beta'};
            if 1-abs(state'*model.transition.(names{branch}))^2>second.overlap_tolerance
                error('Transition overlap is below the explicit tolerance.');
            end
            correction=zeros(cluster_system.bath_dimension);
            for a=1:numel(fields)
                for b=1:numel(fields)
                    weight=0;
                    for n=outside(:).'
                        weight=weight+(state'*ops{a}*vectors(:,n))*...
                                      (vectors(:,n)'*ops{b}*state)/gaps(n);
                    end
                    correction=correction+weight*fields{a}*fields{b};
                end
            end
            if branch==1
                left=left+correction;
            else
                right=right+correction;
            end
        end
        left=(left+left')/2;
        right=(right+right')/2;

    end
    H=struct('kind','cce_projected','left',left,'right',right,...
             'formalism',cluster_system.formalism,'physical_ids',ids,...
             'dissipator',[],'swap_matrix',[]);
    H.hamiltonian_abs_error_bound=meta.roundoff_bound;
    H.hamiltonian_norm=norm(ham,'fro');
    if ~isempty(jumps)||strcmp(cluster_system.formalism,'zeeman-liouv')
        [H.dissipator,jump_records]=dissipator(cluster_system,jumps,va,vb);

        % Map each column-major bath entry to its transposed position
        labels=reshape(1:db^2,db,db).';
        H.swap_matrix=sparse(1:db^2,labels(:),1,db^2,db^2);

    end
    R=[];
    K=[];
elseif strcmp(cluster_system.formalism,'zeeman-liouv')

    % Column-major vectorisation fixes the transposes and the relaxation sign
    memory_guard(parameters,16*((12+3*numel(ids))*d^2+8*d^4));
    H=kron(eye(d),ham)-kron(ham.',eye(d));
    [R,jump_records]=dissipator(cluster_system,jumps,[],[]);
    K=zeros(size(H));
else
    H=ham;
    R=zeros(size(H));
    K=zeros(size(H));

end

% Retain the physical master-generator evidence with its completed raw task
if exist('jump_records','var')
    ha=ham;
    hb=ham;
    if projected
        ha=H.left;
        hb=H.right;
    end
    parameters.cce_aux.generator=struct('physical_ids',ids,'dimension',d,...
        'bath_dimension',db,'observed_dimension',cluster_system.observed_dimension,...
        'projected',projected,'h_alpha',ha,'h_beta',hb,'jumps',jump_records,...
        'vectorization','column','evolution','exp(-1i*L*t)',...
        'rate_units','s^-1','hamiltonian',meta);
end

% Cache only exact native inputs and the corresponding sequence auxiliaries
hits=0;
misses=0;
if ~isempty(cache)
    hits=cache.hits;
    misses=cache.misses;
end
cache=struct('signature',signature,'H',H,'R',R,'K',K,...
             'parameters',parameters,'hits',hits,'misses',misses+1);

end

% Construct the same retained Hamiltonian without invoking the old API
function [H,meta]=matrix(system,parameters)

% Resolve only the construction backend; physical tensors are already in rad/s
model=system.model;
spin_config=system.configuration;
bath_dimension=system.bath_dimension;
backend=parameters.cce.backend;
if strcmp(backend,'auto')
    if model.direct_supported
        backend='direct';
    else
        backend='stock';
    end
end
if strcmp(backend,'direct')&&~model.direct_supported
    error('cce_generators:backend','The requested direct approximation is unsupported.');
end
ids=system.physical_ids;
dims=system.multiplicities;
dimension=system.dimension;
memory_guard(parameters,16*dimension^2*(12+3*numel(ids)));

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
        tensor=model.zeeman{id};
        strength=model.spin_system.inter.zeeman.strength{id};
        field=zeros(1,3);
        isotropic=trace(model.spin_system.inter.zeeman.matrix{id})/3;
        anisotropic=tensor-isotropic*eye(3);
        [~,rank_one,rank_two]=mat2sphten(model.spin_system.inter.zeeman.matrix{id});
        if norm(rank_one,2)<=model.spin_system.tols.liouv_zero&&...
           norm(rank_two,2)<=model.spin_system.tols.liouv_zero
            anisotropic=zeros(3);
        end
        switch strength
            case 'full'
                if abs(isotropic)<=model.spin_system.tols.liouv_zero
                    isotropic=0;
                end
                field=anisotropic(3,:)+[0 0 isotropic];
            case 'secular'
                isotropic=isotropic-model.spin_system.inter.basefrqs(id);
                if abs(isotropic)<=model.spin_system.tols.liouv_zero
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
                term=tensor(axis,other)*(ops{n,axis}*ops{n,other});
                H=H+term;
                sum_norm=sum_norm+norm(term,'fro');
                operations=operations+dimension;
            end
        end

        % Sum each physical pair once, using combined ordered-cell tensors
        for k=n+1:numel(ids)
            tensor=model.retained{id,ids(k)};
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
    system=model.spin_system;
    if isfield(system,'bas')
        system=rmfield(system,'bas');
    end
    system=kill_spin(system,setdiff(model.physical_ids,ids));
    basis_parameters.formalism='zeeman-hilb';
    basis_parameters.approximation={'none'};
    system=basis(system,basis_parameters);
    system=assume(system,model.assumptions{:});
    [H,Q]=hamiltonian(system);
    H=H+orientation(Q,model.orientation);

    % Reorder the stock ascending physical indices into observed-first order
    sorted_ids=sort(ids);
    sorted_dims=model.multiplicities(sorted_ids);
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
    if ~model.direct_supported&&~isempty(spin_config.exterior_ids)
        error('exterior mean fields require validated pair-retention semantics.');
    end
    for n=1:numel(ids)
        for k=1:numel(spin_config.exterior_ids)
            field=model.retained{ids(n),spin_config.exterior_ids(k)}*...
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
                    dimension*model.canonical_roundoff;
if meta.antihermitian_residual>2*meta.roundoff_bound
    error('Hamiltonian has a non-roundoff anti-Hermitian residual.');
end
H=(H+H')/2;

% Return ordering, operators, backend, and approximation provenance
meta.physical_ids=ids;
meta.dimensions=dims;
meta.dimension=dimension;
meta.observed_dimension=prod(model.multiplicities(model.observed_ids));
meta.bath_dimension=bath_dimension;
meta.operators=ops;
meta.backend=backend;
meta.fallback_reason=model.fallback_reason;
meta.source_hash=model.source_hash;
meta.omitted_pairs=model.omitted_pairs;
meta.observable_budget=model.numerics.atol+model.numerics.rtol*norm(H,'fro');

end

% Assemble loss and jump terms for conventional or conditional density blocks
function [R,records]=dissipator(system,jumps,va,vb)

% Select the density dimension after the caller's allocation preflight
ids=system.physical_ids;
d=system.dimension;
projected=~isempty(va);
if projected
    n=system.bath_dimension;
else
    n=d;
end
R=zeros(n^2);
records=struct('support',{},'rate',{},'left',{},'right',{});
for k=1:numel(jumps)

    % A supported jump is included only when its entire support is retained
    jump=jumps(k);
    if jump.rate==0
        continue;
    end
    included=ismember(jump.support,ids);
    if ~any(included)
        if any(ismember(jump.support,system.configuration.exterior_ids))
            error('cce:exteriorJump','An active exterior jump cannot drive a frozen mean field.');
        end
        continue
    end
    if ~all(included)
        error('cce:jumpSupport','A correlated jump crosses this cluster; group its support.');
    end
    if ischar(jump.operator)
        labels={'L+','L-','Lx','Ly','Lz'};
        fields={'p','m','x','y','z'};
        field=fields{strcmp(jump.operator,labels)};
        local=system.local_operators{find(ids==jump.support,1)}.(field);
    else
        local=jump.operator;
    end

    % Embed the declared jump and inspect conditional-manifold closure
    op=embed_jump(local,jump.support,ids,system.multiplicities);
    if projected
        left=va'*op*va;
        right=vb'*op*vb;
        loss=op'*op;
        leakage=norm(op*va-va*left,'fro')+norm(op*vb-vb*right,'fro');
        loss_leakage=norm(loss*va-va*(va'*loss*va),'fro')+...
                     norm(loss*vb-vb*(vb'*loss*vb),'fro');
        if leakage>100*eps*d*norm(op,'fro')||loss_leakage>100*eps*d*norm(loss,'fro')
            error('cce:capability','A jump must preserve both projected manifolds.');
        end
    else
        left=op;
        right=op;
    end

    % Accumulate the GKSL gain and both anticommutator loss terms
    R=R+jump.rate*(kron(conj(right),left)-...
        0.5*kron(eye(n),left'*left)-0.5*kron((right'*right).',eye(n)));
    records(end+1)=struct('support',jump.support,'rate',jump.rate,...
                          'left',left,'right',right); %#ok<AGROW>

end

end

% Embed an ordered multi-spin jump without assuming consecutive tensor slots
function op=embed_jump(local,support,ids,dims)

% Reorder a support-first tensor product into retained physical-spin order
[~,positions]=ismember(support,ids);
rest=setdiff(1:numel(ids),positions,'stable');
order=[positions rest];
ordered_dims=dims(order);
if isscalar(dims)
    op=local;
    return
end
labels=reshape(1:prod(dims),fliplr(ordered_dims));
[~,inverse]=sort(order);
labels=permute(labels,numel(dims)+1-fliplr(inverse));
labels=labels(:);
source=kron(local,eye(prod(dims(rest))));
op=source(labels,labels);

end

% Embed one local matrix in the declared retained tensor-product order
function op=embed(local,position,dims)

% Form the ordered product with identities on all other spins
op=1;
for k=1:numel(dims)
    if k==position
        op=kron(op,local);
    else
        op=kron(op,eye(dims(k)));
    end
end

end

% Select one local Cartesian spin operator without repeated full embeddings
function op=pauli_axis(dimension,axis)
local=pauli(dimension);
op=local.(axis);
end

% Reject allocations exceeding the declared retained-space budget
function memory_guard(parameters,bytes)
if bytes>parameters.cce.numerics.max_memory_bytes
    error('cce:memory','Generator construction exceeds max_memory_bytes.');
end
end

% Enforce the native ownership and representation contract
function grumble(system,parameters,cache)
if isfield(parameters,'sequence')
    error('cce_generators:sequence','Use finite_drive, refocusing_train, or hahn_echo, not parameters.sequence.');
end
if ~isstruct(system)||~isscalar(system)||~isfield(system,'kind')||...
   ~strcmp(system.kind,'cce_subsystem')||~isstruct(parameters)||~isscalar(parameters)||...
   ~isfield(parameters,'cce')||~all(isfield(parameters.cce,{'approximation','numerics','second_order','backend'}))
    error('cce_generators:input','A native retained subsystem and resolved parameters are required.');
end
if ~ischar(parameters.cce.approximation)||...
   ~ismember(parameters.cce.approximation,{'generalized','projected'})
    error('cce_generators:approximation','Declare generalized or projected evolution.');
end
if ~ismember(system.formalism,{'zeeman-hilb','zeeman-liouv'})
    error('cce_generators:formalism','Only complete Zeeman representations are supported.');
end
if ~isempty(cache)&&(~isstruct(cache)||~isscalar(cache)||...
   ~all(isfield(cache,{'signature','H','R','K','parameters','hits','misses'})))
    error('cce_generators:cache','Use [] or this generator service cache.');
end
if ~ischar(parameters.cce.backend)||~isrow(parameters.cce.backend)||...
   ~ismember(parameters.cce.backend,{'auto','direct','stock'})
    error('cce_generators:backend','Backend must be auto, direct or stock.');
end
if isfield(parameters,'bath_jumps')
    jumps=parameters.bath_jumps;
    if ~isstruct(jumps)
        error('cce_options:jumps','bath_jumps must be a structure array.');
    end
    model=system.model;
    selected=[model.observed_ids model.bath_ids];
    for k=1:numel(jumps)
        jump=jumps(k);
        required={'support','operator','rate'};
        unknown=setdiff(fieldnames(jump),required);
        if ~isempty(unknown)
            error('cce_options:unknown','Unknown bath_jumps field: %s.',unknown{1});
        end
        if ~all(isfield(jump,required))
            error('cce_options:jumps','Every bath jump needs support, operator and rate.');
        end
        ids=jump.support;
        if ~isnumeric(ids)||~isreal(ids)||~isrow(ids)||isempty(ids)||...
           any(~isfinite(ids))||any(ids~=fix(ids))||any(ids<1)||...
           any(ids>model.spin_system.comp.nspins)||numel(unique(ids))~=numel(ids)
            error('cce_options:spins','Invalid physical spin-ID row in bath_jumps.support.');
        end
        if ~all(ismember(ids,selected))
            error('cce_options:jumps','Jump support must belong to selected observed or bath spins.');
        end
        if ~isnumeric(jump.rate)||~isreal(jump.rate)||~isscalar(jump.rate)||...
           ~isfinite(jump.rate)||jump.rate<0
            error('cce_options:value','Invalid nonnegative value for bath_jumps.rate.');
        end
        if ischar(jump.operator)
            if ~isrow(jump.operator)||~ismember(jump.operator,{'L+','L-','Lx','Ly','Lz'})
                error('cce_options:choice','Unsupported bath_jumps.operator.');
            end
            if numel(ids)~=1
                error('cce_options:jumps','A named local jump needs one supporting spin.');
            end
        elseif ~isnumeric(jump.operator)||any(~isfinite(jump.operator),'all')||...
               ~isequal(size(jump.operator),[prod(model.multiplicities(ids)) prod(model.multiplicities(ids))])
            error('cce_options:jumps','A jump matrix must match its ordered support.');
        end
    end
    if ~isempty(jumps)&&~strcmp(parameters.cce.approximation,'projected')&&...
       ~strcmp(system.formalism,'zeeman-liouv')
        error('cce_generators:formalism','Dissipation requires zeeman-liouv.');
    end
end
end


