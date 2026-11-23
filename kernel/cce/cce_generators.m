% Construct retained-cluster generators in Spinach angular-frequency units.
% Syntax:
%
%   [H,R,K,spin_system,cache]=cce_generators(spin_system,parameters,cache)
%
% spin_system is cce_subsystem output; parameters is a resolved native
% experiment structure. Conventional generators obey L=H+1i*R+1i*K and
% exp(-1i*L*t). R is the GKSL dissipator in s^-1; K is zero. Dissipation
% requires zeeman-liouv, except for projected branches. A projected
% dissipator acts on an alpha/beta coherence block, not a density. Its
% swapped branch is S*conj(R)*S', where S transposes the bath indices.
% parameters.bath_jumps uses support, operator, and rate; supported spin
% descriptors are L+, L-, Lx, Ly, and Lz, or explicit ordered matrices.
% Zero-rate jumps are omitted before closure checks; active split supports
% and exterior jumps driving a frozen mean field are refused, except that
% periodic_segments explicitly freezes the exterior and omits exterior jumps.
% Projected H is the cell pair {left,right} of branch Hamiltonians, R is
% the dissipator acting on column-vectorised bath coherence, and K is [].
% cce_evolution takes them as {left,right,R}. A refocusing pulse swaps
% the two branches and reverses the projected dissipator.
% Conditional midpoint Hahn projection omits the isolated register terms:
% their two scalar phases cancel exactly between equal swapped intervals.
% parameters.cce.backend selects auto/direct/stock. The returned spin
% system carries the generator data in its cce subfield: roundoff_bound,
% backend, and operators, with projection_leakage and hamiltonian_norm for
% projected generators. Master-equation cce.generator retains the physical
% Hamiltonians, spin operators, and included jump pairs for completed
% raw-task evidence.
% cache is [] or the separate exact generator cache returned here.
%
% Parameters:
%
%   spin_system - retained subsystem returned by cce_subsystem
%   parameters  - resolved native experiment and numerical options
%   cache       - [] or the exact generator cache returned previously
%
% Outputs:
%
%   H,R,K       - retained generators, or the projected branch pair in H
%   spin_system - retained subsystem with the generator data in cce
%   cache       - exact-input generator cache with hit and miss counts
%
% ilya.kuprov@weizmann.ac.il

function [H,R,K,spin_system,cache]=cce_generators(spin_system,parameters,cache)

% Validate ownership before allocating retained-space matrices
grumble(spin_system,parameters,cache);
ids=spin_system.cce.ids;
d=spin_system.cce.dimension;
db=spin_system.cce.bath_dimension;
projected=strcmp(parameters.cce.approximation,'projected');
periodic=strcmp(parameters.cce.expansion,'periodic_segments');
jumps=[];
if isfield(parameters,'bath_jumps')
    jumps=parameters.bath_jumps;
end

% Exact input comparison prevents stale orientation or configuration reuse
stale=intersect(fieldnames(spin_system.cce),{'roundoff_bound','backend','operators',...
                'projection_leakage','hamiltonian_norm','generator'});
spin_system.cce=rmfield(spin_system.cce,stale);
signature=struct('spin_system',spin_system,'parameters',parameters);
if ~isempty(cache)&&isequaln(cache.signature,signature)
    H=cache.H;
    R=cache.R;
    K=cache.K;
    spin_system=cache.spin_system;
    cache.hits=cache.hits+1;
    return
end

% Count resident jump operators before Hamiltonian or dissipator allocation
if ~isempty(jumps)||(projected&&strcmp(spin_system.bas.formalism,'zeeman-liouv'))
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
    memory_guard(spin_system,16*((12+3*numel(ids))*d^2+8*n^4+2*included*d^2));
end

% Preserve the canonical direct and independent stock construction paths
[ham,meta]=matrix(spin_system,parameters);
if projected

    % Conditional branches use the declared observed transition, never a guess
    transition=spin_system.cce.transition;
    va=kron(transition.alpha,eye(db));
    vb=kron(transition.beta,eye(db));
    left=va'*ham*va;
    right=vb'*ham*vb;
    leakage=max(norm(ham*va-va*left,'fro'),norm(ham*vb-vb*right,'fro'));
    if ~isempty(jumps)&&leakage>d*eps*norm(ham,'fro')
        error('cce:population_transfer','Projected master dynamics requires invariant observed manifolds.');
    end
    second=parameters.cce.second_order.enabled;
    if second

        % Resolve isolated observed eigenvectors against the compiled transition
        configuration=rmfield(spin_system.cce.configuration,{'exterior_ids','exterior_expectations'});
        [observed_system,~]=cce_subsystem(spin_system,parameters,[],configuration,[]);
        observed_system.cce.configuration.exterior_ids=[];
        observed_system.cce.configuration.exterior_expectations=zeros(0,3);
        [observed_hamiltonian,~]=matrix(observed_system,parameters);
        [vectors,energies]=eig(full((observed_hamiltonian+observed_hamiltonian')/2),'vector');
        [energies,order]=sort(real(energies));
        vectors=vectors(:,order);
        [~,selected(1)]=max(abs(vectors'*spin_system.cce.transition.alpha).^2);
        [~,selected(2)]=max(abs(vectors'*spin_system.cce.transition.beta).^2);
        second=parameters.cce.second_order;
        observed=parameters.observed;
        ops=cell(1,3*numel(observed));
        fields=cell(size(ops));
        axes={'x','y','z'};

        % Resolve both directions of each retained pair
        [pair_index,tensors]=cce_pair_index(spin_system);

        for c=1:numel(observed)
            for a=1:3
                slot=3*(c-1)+a;
                ops{slot}=embed(pauli_axis(spin_system.comp.mults(observed(c)),axes{a}),...
                               c,spin_system.comp.mults(observed));
                field=zeros(spin_system.cce.bath_dimension);
                for b=1:numel(spin_system.cce.sites)
                    tensor=tensors{1+pair_index(observed(c),spin_system.cce.sites(b))};
                    for axis=1:3
                        field=field+tensor(a,axis)*embed(...
                            pauli_axis(spin_system.comp.mults(spin_system.cce.ids(numel(observed)+b)),axes{axis}),...
                            b,spin_system.comp.mults(spin_system.cce.ids(numel(observed)+1:end)));
                    end
                end
                for b=1:numel(spin_system.cce.configuration.exterior_ids)
                    field=field+eye(spin_system.cce.bath_dimension)*(...
                        tensors{1+pair_index(observed(c),spin_system.cce.configuration.exterior_ids(b))}(a,:)*...
                        spin_system.cce.configuration.exterior_expectations(b,:).');
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
            if 1-abs(state'*spin_system.cce.transition.(names{branch}))^2>second.overlap_tolerance
                error('Transition overlap is below the explicit tolerance.');
            end
            correction=zeros(spin_system.cce.bath_dimension);
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
    H={left,right};
    R=[];
    K=[];
    if ~isempty(jumps)||strcmp(spin_system.bas.formalism,'zeeman-liouv')
        [R,jump_records]=dissipator(spin_system,jumps,va,vb);
        if periodic
            H={sparse(left),sparse(right)};
            R=sparse(R);
        end
    end
elseif strcmp(spin_system.bas.formalism,'zeeman-liouv')

    % Column-major vectorisation fixes the transposes and the relaxation sign
    memory_guard(spin_system,16*((12+3*numel(ids))*d^2+8*d^4));
    H=kron(eye(d),ham)-kron(ham.',eye(d));
    [R,jump_records]=dissipator(spin_system,jumps,[],[]);
    K=zeros(size(H));
else
    H=ham;
    R=zeros(size(H));
    K=zeros(size(H));

end

% Store the generator data in the spin system for the pulse sequence
spin_system.cce.roundoff_bound=meta.roundoff_bound;
spin_system.cce.backend=meta.backend;
spin_system.cce.operators=meta.operators;
if projected
    spin_system.cce.projection_leakage=leakage;
    spin_system.cce.hamiltonian_norm=norm(ham,'fro');
end

% Retain the physical master-generator evidence with its completed raw task
if exist('jump_records','var')
    ha=ham;
    hb=ham;
    if projected
        ha=H{1};
        hb=H{2};
    end
    spin_system.cce.generator=struct('physical_ids',ids,'dimension',d,...
        'bath_dimension',db,'observed_dimension',spin_system.cce.observed_dimension,...
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
             'spin_system',spin_system,'hits',hits,'misses',misses+1);

end

% Construct the same retained Hamiltonian without invoking the old API
function [H,meta]=matrix(system,parameters)

% Resolve only the construction backend; physical tensors are already in rad/s
spin_config=system.cce.configuration;
bath_dimension=system.cce.bath_dimension;
backend=parameters.cce.backend;
if strcmp(backend,'auto')
    if system.cce.direct_supported
        backend='direct';
    else
        backend='stock';
    end
end
ids=system.cce.ids;
conditional=isfield(parameters.cce,'projection')&&strcmp(parameters.cce.projection,'conditional');
observed=ismember(ids,parameters.observed);
locals=~(conditional&observed);
dims=system.comp.mults(system.cce.ids);
dimension=system.cce.dimension;
memory_guard(system,16*dimension^2*(12+3*numel(ids)));

% Resolve both directions of each retained pair
[pair_index,tensors]=cce_pair_index(system);

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

        % The isolated register phase cancels in the symmetric Hahn sequence
        id=ids(n);
        if locals(n)

            % Select Zeeman components and rotating-frame carrier convention
            tensor=system.cce.zeeman{id};
            strength=system.inter.zeeman.strength{id};
            field=zeros(1,3);
            isotropic=trace(system.inter.zeeman.matrix{id})/3;
            anisotropic=tensor-isotropic*eye(3);
            [~,rank_one,rank_two]=mat2sphten(system.inter.zeeman.matrix{id});
            if norm(rank_one,2)<=system.tols.liouv_zero&&...
               norm(rank_two,2)<=system.tols.liouv_zero
                anisotropic=zeros(3);
            end
            switch strength
                case 'full'
                    if abs(isotropic)<=system.tols.liouv_zero
                        isotropic=0;
                    end
                    field=anisotropic(3,:)+[0 0 isotropic];
                case 'secular'
                    isotropic=isotropic-system.inter.basefrqs(id);
                    if abs(isotropic)<=system.tols.liouv_zero
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
            tensor=system.cce.quadratic{id};
            strength=system.cce.quadratic_strength{id};
            switch strength
                case 'strong'

                    % Keep the full symmetric traceless tensor
                case 'secular'
                    tensor=diag([-tensor(3,3)/2 -tensor(3,3)/2 tensor(3,3)]);
                case 'zz'
                    if norm(tensor,2)>2*pi*system.tols.inter_cutoff
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

        end

        % Sum each physical pair once, using combined ordered-cell tensors
        for k=n+1:numel(ids)
            if conditional&&observed(n)&&observed(k)
                continue
            end
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
    reduced=system;
    if isfield(reduced,'bas')
        reduced=rmfield(reduced,'bas');
    end
    reduced=kill_spin(reduced,setdiff(1:system.comp.nspins,ids));
    if conditional

        % Remove only pure-register terms before stock Hamiltonian construction
        register=find(ismember(sort(ids),parameters.observed));
        reduced.inter.zeeman.matrix(register)={zeros(3)};
        reduced.inter.basefrqs(register)=0;
        keys=reshape(cell2mat(reduced.inter.coupling.matrix(:,1:2)),[],2);
        reduced.inter.coupling.matrix(all(ismember(keys,register),2),:)=[];
        reduced.inter.giant.coeff(register)={[]};

    end
    basis_parameters.formalism='zeeman-hilb';
    basis_parameters.approximation={'none'};
    reduced=basis(reduced,basis_parameters);
    reduced=assume(reduced,system.cce.assumptions{:});
    [H,Q]=hamiltonian(reduced);
    H=H+orientation(Q,parameters.orientation);

    % Reorder the stock ascending physical indices into observed-first order
    sorted_ids=sort(ids);
    sorted_dims=system.comp.mults(sorted_ids);
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
canonical_roundoff=system.cce.canonical_roundoff;
if conditional

    % Only retained bath Zeeman terms enter the centred Hamiltonian
    canonical_roundoff=0;
    for id=ids(locals)
        canonical_roundoff=canonical_roundoff+27*eps*...
            (norm(system.inter.zeeman.matrix{id},'fro')+...
             abs(system.inter.basefrqs(id)));
    end

end
meta.roundoff_bound=roundoff*max(sum_norm,realmin)+dimension*canonical_roundoff;
if meta.antihermitian_residual>2*meta.roundoff_bound
    error('Hamiltonian has a non-roundoff anti-Hermitian residual.');
end
H=(H+H')/2;

% Return ordering, operators, backend, and approximation provenance
meta.physical_ids=ids;
meta.dimensions=dims;
meta.dimension=dimension;
meta.observed_dimension=prod(system.comp.mults(parameters.observed));
meta.bath_dimension=bath_dimension;
meta.operators=ops;
meta.backend=backend;
meta.fallback_reason=system.cce.fallback_reason;
meta.source_hash=system.cce.source_hash;
meta.omitted_pairs=system.cce.omitted_pairs;
meta.observable_budget=system.tols.cce_atol+system.tols.cce_rtol*norm(H,'fro');

end

% Assemble loss and jump terms for conventional or conditional density blocks
function [R,records]=dissipator(system,jumps,va,vb)

% Select the density dimension after the caller's allocation preflight
ids=system.cce.ids;
d=system.cce.dimension;
projected=~isempty(va);
if projected
    n=system.cce.bath_dimension;
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
    if ~all(ismember(jump.support,ids))
        continue
    end
    if ischar(jump.operator)
        labels={'L+','L-','Lx','Ly','Lz'};
        fields={'p','m','x','y','z'};
        field=fields{strcmp(jump.operator,labels)};
        local=system.cce.local_operators{find(ids==jump.support,1)}.(field);
    else
        local=jump.operator;
    end

    % Embed the declared jump and inspect conditional-manifold closure
    op=embed_jump(local,jump.support,ids,system.comp.mults(system.cce.ids));
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
function memory_guard(spin_system,bytes)
if bytes>spin_system.tols.cce_mem
    error('cce:memory','Generator construction exceeds sys.tols.cce_mem.');
end
end

% Enforce the native ownership and representation contract
function grumble(system,parameters,cache)
if isfield(parameters,'sequence')
    error('cce_generators:sequence','Use finite_drive, refocusing_train, or hahn_echo, not parameters.sequence.');
end
if ~isstruct(system)||~isscalar(system)||~isfield(system,'cce')||...
   ~isfield(system.cce,'sites')||~isstruct(parameters)||~isscalar(parameters)||...
   ~isfield(parameters,'cce')||~all(isfield(parameters.cce,{'approximation','second_order','backend'}))
    error('cce_generators:input','A native retained subsystem and resolved parameters are required.');
end
if ~ischar(parameters.cce.approximation)||...
   ~ismember(parameters.cce.approximation,{'generalized','projected'})
    error('cce_generators:approximation','Declare generalized or projected evolution.');
end
if ~ismember(system.bas.formalism,{'zeeman-hilb','zeeman-liouv'})
    error('cce_generators:formalism','Only complete Zeeman representations are supported.');
end
if ~isempty(cache)&&(~isstruct(cache)||~isscalar(cache)||...
   ~all(isfield(cache,{'signature','H','R','K','spin_system','hits','misses'})))
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
    selected=[parameters.observed parameters.cce.bath_spins];
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
           any(ids>system.comp.nspins)||numel(unique(ids))~=numel(ids)
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
               ~isequal(size(jump.operator),[prod(system.comp.mults(ids)) prod(system.comp.mults(ids))])
            error('cce_options:jumps','A jump matrix must match its ordered support.');
        end
        included=ismember(ids,system.cce.ids);
        if jump.rate~=0&&~any(included)&&...
           any(ismember(ids,system.cce.configuration.exterior_ids))&&...
           ~(strcmp(parameters.cce.approximation,'projected')&&...
             strcmp(parameters.cce.expansion,'periodic_segments'))
            error('cce:exteriorJump','An active exterior jump cannot drive a frozen mean field.');
        end
        if jump.rate~=0&&any(included)&&~all(included)
            error('cce:jumpSupport','A correlated jump crosses this cluster; group its support.');
        end
    end
    if ~isempty(jumps)&&strcmp(parameters.cce.approximation,'projected')&&...
       parameters.cce.second_order.enabled
        error('cce:capability','Dissipative virtual-state corrections are unsupported.');
    end
    if ~isempty(jumps)&&~strcmp(parameters.cce.approximation,'projected')&&...
       ~strcmp(system.bas.formalism,'zeeman-liouv')
        error('cce_generators:formalism','Dissipation requires zeeman-liouv.');
    end
end
if strcmp(parameters.cce.approximation,'projected')&&~isfield(system.cce,'transition')
    error('cce_generators:transition','Projected evolution needs a transition.');
end
if strcmp(parameters.cce.backend,'direct')&&~system.cce.direct_supported
    error('cce_generators:backend','The requested direct approximation is unsupported.');
end
if ~system.cce.direct_supported&&isfield(system.cce.configuration,'exterior_ids')&&...
   ~isempty(system.cce.configuration.exterior_ids)
    error('exterior mean fields require validated pair-retention semantics.');
end
end


