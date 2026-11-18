% Ideal Hahn echo; one pi pulse at half the elapsed time.
% Syntax:
%
%           [E,aux]=hahn_echo(spin_system,parameters,H,R,K)
%
% parameters.times contains total sequence durations in seconds. transition
% selects the observed two-state transition; rho0/coil or preparation may
% override the prepared density and detector. normalized defaults to true
% and divides by the actual detected initial reference. nloops is required
% for CPMG. bath_pumped_hahn requires pump.center, halfwidth, units (Hz or
% rad/s), fraction and allow_empty; optional sites/species restrict the bath.
% Ideal refocusing acts on the selected transition about x. H,R,K are the
% conventional generators or explicitly tagged projected branch generators.
%
% output=density returns column-vectorised generalized observed density.
% It requires a physical initial density and does not normalize a receiver.
%
% cce.numerics.propagator='expm' selects direct segment exponentials for
% generalized Hilbert evolution, retaining the segment arithmetic bound.
% The default reuses a Hermitian eigensystem across the time grid.
% Optional refocusing is a scalar descriptor with spins, operator (a
% Hermitian integrated rotation generator), and either fraction in [0,1]
% or time in seconds. Pulses beyond a sample endpoint are omitted. It
% replaces the default observed pi pulse; omitted bath sites are inactive
% only without a static exterior field. Joint observed/bath pulses are refused.
% Projected pulses must preserve or interchange the selected manifolds.
% Alternatively, refocusing.axis and angle specify a physical Cartesian
% rotation. A one-site exterior Cartesian pulse updates its product density
% and subsequent static mean field; generalized Hilbert propagation then uses
% direct segment exponentials. Exterior explicit operator pulses remain refused.
% aux.pulse_time and aux.pulse_applied record sample-specific control timing.
% Optional aux contains spectral/error diagnostics for coherent projected echoes.
%
% Generalized unitary Zeeman-Hilbert echoes retain observed density matrices,
% spectral residuals and absolute forward-error bounds in aux. The response
% remains time-by-detector, or column-vectorized observed density by time.
%
% ilya.kuprov@weizmann.ac.il

function [E,aux]=hahn_echo(spin_system,parameters,H,R,K)

% Preparation and detection are owned by this experiment; projected states are coherence blocks
grumble(spin_system,parameters,H,R,K);
if ~isfield(parameters,'normalized'), parameters.normalized=true; end
density=isfield(parameters,'output')&&strcmp(parameters.output,'density');
if density
    parameters.normalized=false;
    if ~isfield(parameters,'coil'), parameters.coil=eye(spin_system.observed_dimension); end
end
[rho0,coil,reference,meta]=cce_observation(spin_system,parameters);
if ~iscell(coil), coil={coil}; end
projected=isstruct(H); E=complex(zeros(numel(parameters.times),numel(reference)));
aux=struct();
if density
    if projected, error('hahn_echo:density','Projected blocks do not define an observed density.'); end
    E=complex(zeros(numel(parameters.times),spin_system.observed_dimension^2));
end
if projected, generator=H; else, generator=H+1i*R+1i*K; end

% Build ideal observed-spin transition pulses without a global bath basis
model=spin_system.model; dc=spin_system.observed_dimension;
db=spin_system.bath_dimension;
if ~isfield(model,'transition')
    error('hahn_echo:transition','An explicit transition is required for ideal refocusing.');
end
a=model.transition.alpha; b=model.transition.beta;
rotation=expm(-1i*pi*kron((a*b'+b*a')/2,eye(db)));
positions=parameters.times/2; H_after=H;
if isfield(parameters,'refocusing')
    pulse=parameters.refocusing;
    if isfield(pulse,'fraction'), positions=pulse.fraction*parameters.times;
    else, positions=repmat(pulse.time,size(parameters.times)); end
    if isfield(pulse,'axis')
        if strcmp(parameters.cce.meanfield,'static')&&any(ismember(pulse.spins,model.bath_ids))
            parameters.cce.numerics.propagator='expm';
        end
        axis=pulse.axis; if ischar(axis), axis=double('xyz'==axis); end
        generator=zeros(dc*db); local=spin_system; local.formalism='zeeman-hilb';
        for id=intersect(pulse.spins,spin_system.physical_ids,'stable')
            for k=1:3
                labels={'Lx','Ly','Lz'};
                generator=generator+pulse.angle*axis(k)*cce_operator(local,labels{k},id);
            end
        end
        exterior=setdiff(pulse.spins,spin_system.physical_ids);
        if ~isempty(exterior)&&strcmp(parameters.cce.meanfield,'static')
            if projected||~strcmp(spin_system.formalism,'zeeman-hilb')
                error('hahn_echo:exterior','Exterior rotations require generalized Hilbert dynamics.');
            end
            configuration=spin_system.configuration;
            configuration=rmfield(configuration,{'expectations','exterior_ids','exterior_expectations'});
            for id=exterior
                index=find(configuration.bath_ids==id,1); operators=pauli(model.multiplicities(id));
                rotation=expm(-1i*pulse.angle*(axis(1)*operators.x+axis(2)*operators.y+axis(3)*operators.z));
                configuration.bath_rho{index}=rotation*configuration.bath_rho{index}*rotation';
            end
            [after,~]=cce_subsystem(model,spin_system.sites,configuration,[]);
            [H_after,~,~]=cce_generators(after,parameters,[]);
            parameters.cce.numerics.propagator='expm';
        end
    elseif isequal(pulse.spins,model.observed_ids)
        generator=kron(pulse.operator,eye(db));
    else
        position=find(spin_system.physical_ids==pulse.spins,1);
        if isempty(position)
            if strcmp(parameters.cce.meanfield,'static')
                error('hahn_echo:exterior','Exterior pulses require time-dependent mean fields.');
            end
            generator=zeros(dc*db);
        else
            dims=spin_system.multiplicities;
            generator=kron(eye(prod(dims(1:position-1))),...
                           kron(pulse.operator,eye(prod(dims(position+1:end)))));
        end
    end
    rotation=expm(-1i*generator);
end

% Retain spectral manifold propagation for a supported ideal pulse
if projected&&~density&&isempty(H.dissipator)&&...
   ~strcmp(spin_system.formalism,'zeeman-liouv')&&~isfield(parameters,'pump')&&...
   all(isfield(H,{'hamiltonian_abs_error_bound','hamiltonian_norm'}))
    [E,aux]=projected_echo(spin_system,parameters,H,reference,meta,rotation,positions);
    aux.pulse_time=positions; aux.pulse_applied=parameters.times>0&positions<=parameters.times; return
end

% Keep ideal master-equation propagation on the native stock action path
if strcmp(spin_system.formalism,'zeeman-liouv')||(projected&&~isempty(H.dissipator))
    [E,aux]=master_echo(spin_system,parameters,H,R,K,rho0,coil,meta,rotation,positions,density);
    aux.pulse_time=positions; aux.pulse_applied=parameters.times>0&positions<=parameters.times; return
end
fractions=0.5;

% Reuse one Hermitian spectrum and retain its forward-error diagnostics
if ~projected&&strcmp(spin_system.formalism,'zeeman-hilb')&&...
   ~any(R(:))&&~any(K(:))&&~isfield(parameters,'pump')&&...
   isfield(parameters,'cce_aux')
    [E,aux]=generalized_echo(spin_system,parameters,H,H_after,rho0,coil,meta,rotation,positions,density);
    aux.pulse_time=positions; aux.pulse_applied=parameters.times>0&positions<=parameters.times; return
end

% At zero elapsed time no control is applied, matching the reference signal
for n=1:numel(parameters.times)
    rho=rho0; previous=0; swapped=false;
    if parameters.times(n)>0
        for j=1:numel(fractions)
            segment=generator;
            if projected&&swapped
                segment.left=H.right; segment.right=H.left;
                if ~isempty(H.dissipator)
                    segment.dissipator=H.swap_matrix*conj(H.dissipator)*H.swap_matrix';
                end
            end
            duration=(fractions(j)-previous)*parameters.times(n);
            rho=cce_evolution(spin_system,segment,[],rho,duration,1,'final');
            if projected
                swapped=~swapped;
            else
                rho=reshape(rho,dc*db,dc*db); rho=rotation*rho*rotation';
                if strcmp(spin_system.formalism,'zeeman-liouv'), rho=rho(:); end
            end
            previous=fractions(j);
        end
    end
    segment=generator;
    if projected&&swapped
        segment.left=H.right; segment.right=H.left;
        if ~isempty(H.dissipator)
            segment.dissipator=H.swap_matrix*conj(H.dissipator)*H.swap_matrix';
        end
    end
    rho=cce_evolution(spin_system,segment,[],rho,(1-previous)*parameters.times(n),1,'final');
    if density
        rho=reshape(rho,dc*db,dc*db); observed=complex(zeros(dc));
        for a=1:dc
            for b=1:dc
                observed(a,b)=trace(rho((a-1)*db+(1:db),(b-1)*db+(1:db)));
            end
        end
        E(n,:)=reshape(observed,1,[]); continue
    end
    for k=1:numel(coil)
        if isvector(rho), value=coil{k}'*rho;
        else, value=trace(coil{k}'*rho); end
        if meta.normalized, value=value/reference(k); end
        E(n,k)=value;
    end
end

end

% Validate the registered sequence contract before preparation
function grumble(spin_system,parameters,H,R,K)
if ~isstruct(spin_system)||~isfield(spin_system,'kind')||...
   ~strcmp(spin_system.kind,'cce_subsystem')||~isstruct(parameters)||...
   ~isfield(parameters,'times')||~isnumeric(parameters.times)||...
   ~isreal(parameters.times)||~isvector(parameters.times)||...
   any(~isfinite(parameters.times))||any(parameters.times<0)||any(diff(parameters.times)<0)
    error('hahn_echo:input','A native retained subsystem and increasing nonnegative times are required.');
end
if isfield(parameters,'sequence')
    error('hahn_echo:sequence','Use finite_drive, refocusing_train, or this named experiment, not parameters.sequence.');
end
if ~(isnumeric(H)||isstruct(H))||~isnumeric(R)||~isnumeric(K)
    error('hahn_echo:generators','Supply conventional or tagged projected generators.');
end
if isfield(parameters,'refocusing')
    pulse=parameters.refocusing;
    if ~isstruct(pulse)||~isscalar(pulse)||...
       ~isfield(pulse,'spins')||...
       ~xor(isfield(pulse,'operator'),all(isfield(pulse,{'axis','angle'})))||...
       ~xor(isfield(pulse,'fraction'),isfield(pulse,'time'))||...
       ~isempty(setdiff(fieldnames(pulse),{'spins','operator','axis','angle','fraction','time'}))
        error('hahn_echo:refocusing','Supply spins, operator, and either fraction or time for one ideal refocusing pulse.');
    end
    ids=pulse.spins;
    if ~isnumeric(ids)||~isreal(ids)||~isrow(ids)||isempty(ids)||...
       any(~ismember(ids,spin_system.model.physical_ids))||numel(unique(ids))~=numel(ids)||...
       ~(isscalar(ids)||isequal(ids,spin_system.model.observed_ids))
        error('hahn_echo:operator_target','An operator pulse addresses the whole observed system or one physical site.');
    end
    d=prod(spin_system.model.multiplicities(ids));
    if isfield(pulse,'operator')&&(~isnumeric(pulse.operator)||~isequal(size(pulse.operator),[d d])||...
       any(~isfinite(pulse.operator),'all')||...
       norm(pulse.operator-pulse.operator','fro')>d*eps*max(1,norm(pulse.operator,'fro')))
        error('hahn_echo:refocusing','The integrated generator must be Hermitian.');
    end
    if isfield(pulse,'axis')
        axis=pulse.axis;
        if ischar(axis)&&ismember(axis,{'x','y','z'}), axis=double('xyz'==axis); end
        if ~isnumeric(axis)||~isreal(axis)||~isequal(size(axis),[1 3])||...
           any(~isfinite(axis))||abs(norm(axis)-1)>3*eps||...
           ~isnumeric(pulse.angle)||~isreal(pulse.angle)||~isscalar(pulse.angle)||~isfinite(pulse.angle)
            error('hahn_echo:refocusing','Supply a unit Cartesian axis and finite real rotation angle.');
        end
    end
    if isfield(pulse,'fraction'), position=pulse.fraction; upper=1;
    else, position=pulse.time; upper=Inf; end
    if ~isnumeric(position)||~isreal(position)||~isscalar(position)||...
       ~isfinite(position)||position<0||position>upper
        error('hahn_echo:refocusing','Use a fraction in [0,1] or nonnegative time in seconds.');
    end
end
if isfield(parameters.cce.numerics,'propagator')&&...
   strcmp(parameters.cce.numerics.propagator,'expm')&&...
   (isstruct(H)||~strcmp(spin_system.formalism,'zeeman-hilb')||...
    any(R(:))||any(K(:))||isfield(parameters,'pump'))
    error('hahn_echo:propagator','Direct segment exponentials require closed generalized Hilbert dynamics.');
end

end

% Propagate the two invariant branches once per interval; retain spectral bounds
function [E,aux]=projected_echo(system,parameters,H,reference,meta,pulse,positions)

% Resolve the branch spectra once, not once per observation time
[va,ea]=eig(full((H.left+H.left')/2),'vector');
[vb,eb]=eig(full((H.right+H.right')/2),'vector');
ra=norm(H.left*va-va.*ea.','fro'); rb=norm(H.right*vb-vb.*eb.','fro');
oa=norm(va'*va-eye(size(H.left)),'fro'); ob=norm(vb'*vb-eye(size(H.right)),'fro');
residual=max(ra,rb); orthogonality=max(oa,ob); db=system.bath_dimension;

% Resolve whether the physical pulse preserves or interchanges the manifolds
transition=system.model.transition;
pa=kron(transition.alpha,eye(db)); pb=kron(transition.beta,eye(db));
aa=pa'*pulse*pa; ab=pa'*pulse*pb; ba=pb'*pulse*pa; bb=pb'*pulse*pb;
tolerance=system.dimension^3*eps;
if norm(aa,'fro')+norm(bb,'fro')<=tolerance
    swap=true; pair={ba,ab};
elseif norm(ab,'fro')+norm(ba,'fro')<=tolerance
    swap=false; pair={aa,bb};
else
    error('hahn_echo:manifolds','The pulse must preserve or interchange the selected spin manifolds.');
end
if max(norm(pair{1}'*pair{1}-eye(db),'fro'),norm(pair{2}'*pair{2}-eye(db),'fro'))>tolerance
    error('hahn_echo:manifolds','The pulse must not leak outside the selected spin manifolds.');
end
rho_b=1;
for id=system.sites
    where=find(system.configuration.bath_ids==id,1);
    rho_b=kron(rho_b,system.configuration.bath_rho{where});
end

% Keep the branch-product order and zero-time convention explicit
E=complex(zeros(numel(parameters.times),1)); bounds=zeros(size(E));
for n=1:numel(parameters.times)
    time=parameters.times(n);
    if time==0
        value=trace(rho_b);
    else
        ua=eye(db); ub=eye(db); previous=0; swapped=false;
        if positions(n)<=time
            dt=positions(n);
            ua=((va.*exp(-1i*dt*ea).')*va')*ua;
            ub=((vb.*exp(-1i*dt*eb).')*vb')*ub;
            ua=pair{1}*ua; ub=pair{2}*ub; swapped=swap; previous=dt;
        end
        dt=time-previous;
        left=(va.*exp(-1i*dt*ea).')*va'; right=(vb.*exp(-1i*dt*eb).')*vb';
        if swapped, ua=right*ua; ub=left*ub;
        else, ua=left*ua; ub=right*ub; end
        value=trace(ua*rho_b*ub');
    end
    E(n)=value;
    operations=20*db^3; gamma=operations*eps/(1-operations*eps);
    bounds(n)=2*time*(residual+H.hamiltonian_abs_error_bound)+8*orthogonality+...
              gamma*(1+time*H.hamiltonian_norm);
end

% A selected-transition detector scales only unnormalised responses
if ~meta.normalized
    E=E*reference; bounds=bounds*abs(reference);
end
aux=struct('abs_error_bound',bounds,'eigen_residual',residual,...
           'orthogonality_residual',orthogonality,...
           'hamiltonian_bound',H.hamiltonian_abs_error_bound,...
           'backend',parameters.cce_aux.backend);

end

% Propagate a generalized Hahn echo with retained observed density and bounds
function [answer,aux]=generalized_echo(system,parameters,H,H_after,rho0,coil,meta,pulse,positions,density)

% Select direct segments or reuse one Hermitian eigensystem
direct=isfield(parameters.cce.numerics,'propagator')&&...
       strcmp(parameters.cce.numerics.propagator,'expm');
if ~direct
    H=(H+H')/2; [vectors,energies]=eig(full(H),'vector'); energies=real(energies);
    residual=norm(H*vectors-vectors.*energies.','fro');
    orthogonality=norm(vectors'*vectors-eye(size(H)),'fro');
end
dc=system.observed_dimension; db=system.bath_dimension; count=numel(parameters.times);
detectors=cellfun(@(x)x(1:db:end,1:db:end),coil,'UniformOutput',false);
answer=complex(zeros(count,numel(meta.normalization_divisor)));
if density, answer=complex(zeros(count,dc^2)); end
aux.density=complex(zeros(dc,dc,count)); bounds=zeros(size(answer));

% Spectral and detector norms are constant over the requested time mesh
operations=20*size(H,1)^3; gamma=operations*eps/(1-operations*eps);
hnorm=norm(H,'fro'); hbound=parameters.cce_aux.roundoff_bound;
detector_norms=cellfun(@(x)norm(x,2),detectors);
if ~isfield(parameters,'coil')&&~isfield(parameters,'detection'), detector_norms=1; end
divisors=meta.normalization_divisor;

% Compose both free intervals and the physical pi rotation in time order
for n=1:count
    time=parameters.times(n); propagator=eye(size(H));
    if direct
        rho=rho0; durations=0; scale=0;
        if time>0&&positions(n)<=time, durations=[0 positions(n) time-positions(n)];
        elseif time>0, durations=[0 time]; end
        for k=1:numel(durations)
            rotation=eye(size(H));
            if k==2&&positions(n)<=time, rotation=pulse; end
            segment=H; if k==3, segment=H_after; end
            propagator=rotation*expm(-1i*durations(k)*segment);
            rho=propagator*rho*propagator';
            scale=scale+1+abs(durations(k))*norm(full(segment),2);
        end
    else
        if time>0
            previous=0;
            if positions(n)<=time
                free=(vectors.*exp(-1i*energies*positions(n)).')*vectors';
                propagator=free*propagator; propagator=pulse*propagator; previous=positions(n);
            end
            free=(vectors.*exp(-1i*energies*(time-previous)).')*vectors';
            propagator=free*propagator;
        end
        rho=propagator*rho0*propagator';
    end
    observed=complex(zeros(dc));
    for k=1:db
        indices=k:db:dc*db; observed=observed+rho(indices,indices);
    end
    aux.density(:,:,n)=observed;
    if density
        answer(n,:)=observed(:).';
    else
        if ~isfield(parameters,'coil')&&~isfield(parameters,'detection')
            answer(n)=exp(-1i*meta.receiver_phase)*(meta.transition.alpha'*...
                      observed*meta.transition.beta)/divisors;
        else
            for k=1:numel(detectors), answer(n,k)=trace(detectors{k}'*observed)/divisors(k); end
        end
    end

    % Bound spectral residuals, tensor assembly and floating-point products
    if direct
        bound=16*size(H,1)^3*eps*scale;
    else
        bound=2*time*(residual+hbound)+...
              8*orthogonality+gamma*(1+time*hnorm);
    end
    if density
        bounds(n,:)=bound;
    else
        for k=1:numel(meta.normalization_divisor)
            bounds(n,k)=bound*detector_norms(k)/abs(divisors(k));
        end
    end
end

% Keep generator and propagation evidence alongside the natural response
aux.abs_error_bound=bounds;
if ~direct
    aux.eigen_residual=residual; aux.orthogonality_residual=orthogonality;
    aux.hamiltonian_bound=parameters.cce_aux.roundoff_bound;
    aux.backend=parameters.cce_aux.backend;
end

end

% Propagate an ideal pulse between native master-equation intervals
function [answer,aux]=master_echo(system,p,H,R,K,rho0,coils,meta,pulse,positions,density)

% Preserve column-major density ordering and the Spinach generator sign
dc=system.observed_dimension; db=system.bath_dimension; projected=isstruct(H);
if projected
    G=-1i*(kron(eye(db),H.left)-kron(H.right.',eye(db)));
    if ~isempty(H.dissipator), G=G+H.dissipator; end
    flip=reshape(1:db^2,db,db).'; flip=flip(:); alternate=conj(G(flip,flip));
    a=kron(system.model.transition.alpha,eye(db)); b=kron(system.model.transition.beta,eye(db));
    aa=a'*pulse*a; ab=a'*pulse*b; ba=b'*pulse*a; bb=b'*pulse*b;
    tolerance=100*eps*system.dimension;
    if norm(ab,'fro')+norm(ba,'fro')<tolerance
        pair={aa,bb}; swapped=false;
    elseif norm(aa,'fro')+norm(bb,'fro')<tolerance
        pair={ba,ab}; swapped=true;
    else
        error('hahn_echo:manifolds','The pulse must preserve or interchange the selected spin manifolds.');
    end
    if norm(pair{1}'*pair{1}-eye(db),'fro')+...
       norm(pair{2}'*pair{2}-eye(db),'fro')>tolerance
        error('hahn_echo:manifolds','The pulse must not leak outside the selected spin manifolds.');
    end
    rho0=1;
    for state=system.bath_rho, rho0=kron(rho0,state{1}); end
    rho0=rho0(:); answer=complex(zeros(numel(p.times),1));
else
    G=-1i*H+R+K;
    if ~isfield(p,'bath_jumps')||isempty(p.bath_jumps), G=sparse(G); end
    answer=complex(zeros(numel(p.times),numel(meta.normalization_divisor)));
    if density, answer=complex(zeros(numel(p.times),dc^2)); end
end
engine=system.model.spin_system; engine.bas.formalism='sphten-liouv'; engine.tols.prop_chop=0;
aux.density=complex(zeros(dc,dc,numel(p.times)));
if projected, aux.density=[]; end

% Use the same stock exponential action on both sides of the ideal rotation
for n=1:numel(p.times)
    time=p.times(n); vector=rho0; generator=G; previous=0;
    if time>0&&positions(n)<=time
        previous=positions(n);
        if isfield(p.cce.numerics,'propagator')&&strcmp(p.cce.numerics.propagator,'stock')
            vector=propagator(engine,1i*G,previous)*vector;
        else
            vector=step(engine,1i*G,vector,previous);
        end
        if projected
            rho=reshape(vector,db,db); rho=pair{1}*rho*pair{2}';
            if swapped, generator=alternate; end
        else
            rho=reshape(vector,dc*db,dc*db); rho=pulse*rho*pulse';
        end
        vector=rho(:);
    end
    if isfield(p.cce.numerics,'propagator')&&strcmp(p.cce.numerics.propagator,'stock')
        vector=propagator(engine,1i*generator,time-previous)*vector;
    else
        vector=step(engine,1i*generator,vector,time-previous);
    end
    if projected
        answer(n)=trace(reshape(vector,db,db));
        if ~meta.normalized, answer(n)=answer(n)*meta.initial_signal; end
    else
        rho=reshape(vector,dc*db,dc*db); observed=complex(zeros(dc));
        for a=1:dc
            for b=1:dc
                observed(a,b)=trace(rho((a-1)*db+(1:db),(b-1)*db+(1:db)));
            end
        end
        aux.density(:,:,n)=observed;
        if density
            answer(n,:)=observed(:).';
        elseif ~isfield(p,'coil')&&~isfield(p,'detection')
            answer(n)=exp(-1i*meta.receiver_phase)*(meta.transition.alpha'*observed*...
                      meta.transition.beta)/meta.normalization_divisor;
        else
            for k=1:numel(coils)
                detector=reshape(coils{k},dc*db,dc*db); detector=detector(1:db:end,1:db:end);
                answer(n,k)=trace(detector'*observed)/meta.normalization_divisor(k);
            end
        end
    end
end

% Retain the master-equation arithmetic estimate, not a physical truncation bound
scale=ones(1,size(answer,2));
if ~projected&&~density
    if isfield(p,'coil')
        for k=1:numel(coils)
            detector=reshape(coils{k},dc*db,dc*db); scale(k)=norm(detector(1:db:end,1:db:end),2);
        end
    end
    scale=scale./abs(meta.normalization_divisor);
elseif projected&&~meta.normalized
    scale=abs(meta.initial_signal);
end
aux.abs_error_bound=(64*eps*size(G,1)*(1+norm(G,1)*p.times(:)))*scale;
aux.error_bound_kind='floating_point_estimate';

end


