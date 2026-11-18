% Internal acquisition for native free evolution and driven experiments.
% Syntax:
%
%     [answer,aux]=cce_acquire(spin_system,parameters,H,R,K)
%
% parameters.times gives elapsed times in seconds. rho0 and coil are explicit
% observed-spin matrices, or transition selects preparation and detection.
% preparation may contain state and detector descriptors. normalized defaults
% to true for signals. output is signal, observable or density; density output
% is the column-vectorised observed density, one row per elapsed time.
% Optional aux retains coherent spectral residuals and absolute error bounds.
% Finite Hilbert controls use segmentwise expm and preserve the original
% 16*d^3*eps*sum(1+abs(dt)*norm(H,2)) arithmetic bound, including the initial
% zero-duration interval.
% Native drive, refocusing, and train declarations are compiled internally.
% Caller-supplied sequence fields are rejected, including empty arrays.
% Drift and relaxation remain active during drives; overlapping drives add.
% Ideal and static-exterior Hilbert segments use cce_segments; pure finite
% drives retain streamed propagation. Projected drives are unsupported.
% H,R,K are conventional angular-frequency generators or explicitly tagged
% projected branches from cce_generators. No pulse or receiver is implicit.
%
% GKSL relaxation uses stock exponential action and its floating-point
% estimate 64*eps*n*(1+norm(L,1)*time), scaled by the detector/reference.
% Master-equation generators retain dense arithmetic. L=H+1i*R+1i*K.
% Projected propagation carries a
% coherence block, not a physical density; density requires generalized CCE.
%
% ilya.kuprov@weizmann.ac.il

function [answer,aux]=cce_acquire(spin_system,parameters,H,R,K)

% Preparation and detection belong to this experiment, not to its context
grumble(spin_system,parameters,H,R,K);
events=cce_pulse_events(parameters);
if ~isfield(parameters,'normalized'), parameters.normalized=true; end
if isfield(parameters,'output')&&strcmp(parameters.output,'density')
    parameters.normalized=false;
    if ~isfield(parameters,'coil'), parameters.coil=eye(spin_system.observed_dimension); end
end
[rho0,coil,reference,meta]=cce_observation(spin_system,parameters);
selected=~isfield(parameters,'coil')&&~isfield(parameters,'detection')&&isfield(meta,'transition');

% Preserve the projected spectral action for a selected-transition FID
if isstruct(H)&&isempty(H.dissipator)&&~strcmp(spin_system.formalism,'zeeman-liouv')&&...
   ~(isfield(parameters,'output')&&strcmp(parameters.output,'density'))&&isempty(events)&&...
   ~isfield(parameters,'coil')&&~isfield(parameters,'detection')&&...
   ~isfield(parameters,'bath_operator')
    [answer,aux]=projected_fid(spin_system,parameters,H,reference,meta); return
end
if ~iscell(coil), detectors={coil}; else, detectors=coil; end
if isstruct(H), generator=H; else, generator=H+1i*R+1i*K; end
density=isfield(parameters,'output')&&strcmp(parameters.output,'density');
if density
    if isstruct(H), error('cce_acquire:density','Projected branches do not return a density matrix.'); end
    answer=complex(zeros(numel(parameters.times),spin_system.observed_dimension^2));
else
    answer=complex(zeros(numel(parameters.times),numel(reference)));
end

% Compile experiment-owned finite drives once for this retained subsystem
segments=struct(); segmented=false;
if ~isempty(events)
    if isstruct(H), error('cce_acquire:projection','Finite controls require generalized dynamics.'); end
    segmented=strcmp(spin_system.formalism,'zeeman-hilb')&&...
              (any(strcmp({events.type},'ideal'))||strcmp(parameters.cce.meanfield,'static')||...
               any(~ismember([events.target],spin_system.physical_ids)));
    if segmented
        segments=cce_segments(spin_system,parameters,H,events);
    else
        [starts,drives]=finite_controls(spin_system,parameters,events);
    end
end

% Reuse the coherent spectrum and retain its measured numerical diagnostics
aux=struct();
controlled=~isempty(events)&&~isstruct(H)&&strcmp(spin_system.formalism,'zeeman-hilb');
if controlled
    aux.abs_error_bound=zeros(size(answer));
    if segmented, aux.segments=segments; end
    aux.density=complex(zeros(spin_system.observed_dimension,spin_system.observed_dimension,numel(parameters.times)));
end
coherent=isempty(events)&&~isstruct(H)&&strcmp(spin_system.formalism,'zeeman-hilb');
projected_me=isstruct(H)&&strcmp(spin_system.formalism,'zeeman-liouv');
liouville=isempty(events)&&(projected_me||(~isstruct(H)&&...
          strcmp(spin_system.formalism,'zeeman-liouv')));
if projected_me
    db=spin_system.bath_dimension;
    generator=kron(eye(db),H.left)-kron(H.right.',eye(db));
    if ~isempty(H.dissipator), generator=generator+1i*H.dissipator; end
    if meta.normalized
        rho0=rho0/meta.initial_signal; reference=1;
    end
end

% Preserve the stock density action and its arithmetic estimate
if liouville
    engine=spin_system.model.spin_system; engine.bas.formalism='sphten-liouv';
    engine.tols.prop_chop=0;
    scale=ones(1,size(answer,2));
    if projected_me
        scale=norm(reshape(full(coil),spin_system.bath_dimension,[]),2);
        if ~meta.normalized, scale=scale*abs(meta.initial_signal); end
    elseif ~density
        if ~selected
            db=spin_system.bath_dimension; dc=spin_system.observed_dimension;
            for k=1:numel(detectors)
                detector=reshape(detectors{k},dc*db,dc*db);
                scale(k)=norm(detector(1:db:end,1:db:end),2);
            end
        end
        if meta.normalized, scale=scale./abs(reference); end
    end
    argument=norm(-1i*generator,1)*parameters.times(:);
    aux.abs_error_bound=(64*eps*size(generator,1)*(1+argument))*scale;
    aux.error_bound_kind='floating_point_estimate';
end
if coherent
    if any(~isfinite(generator),'all')||...
       norm(generator-generator','fro')>size(generator,1)*eps(max(1,norm(generator,'fro')))
        error('cce_evolution:hermiticity','Hilbert generators must be finite and Hermitian.');
    end
    [vectors,energies]=eig(full((generator+generator')/2),'vector');
    aux.eigen_residual=norm(H*vectors-vectors.*energies.','fro');
    aux.orthogonality_residual=norm(vectors'*vectors-eye(size(generator)),'fro');
    aux.hamiltonian_bound=parameters.cce_aux.roundoff_bound;
    aux.backend=parameters.cce_aux.backend;
    operations=12*size(generator,1)^3; gamma=operations*eps/(1-operations*eps);
    bounds=2*parameters.times(:)*(aux.eigen_residual+aux.hamiltonian_bound)+...
           4*aux.orthogonality_residual+...
           gamma*(1+parameters.times(:)*norm(generator,'fro'));
    if density
        scale=ones(1,spin_system.observed_dimension^2);
    else
        scale=cellfun(@(op)norm(op,2),detectors(:).');
        if meta.normalized, scale=scale./abs(reference); end
    end
    aux.abs_error_bound=bounds.*scale;
end

% Every requested elapsed time starts from the same prepared density matrix
for n=1:numel(parameters.times)
    if coherent
        if parameters.times(n)==0
            rho=rho0;
        else
            u=vectors*diag(exp(-1i*energies*parameters.times(n)))*vectors';
            rho=u*rho0*u';
        end
    elseif liouville
        if isfield(parameters.cce.numerics,'propagator')&&...
           strcmp(parameters.cce.numerics.propagator,'stock')
            rho=propagator(engine,generator,parameters.times(n))*rho0;
        else
            action=generator;
            if ~projected_me&&(~isfield(parameters,'bath_jumps')||isempty(parameters.bath_jumps))
                action=sparse(action);
            end
            rho=step(engine,action,rho0,parameters.times(n));
        end
    elseif isempty(events)
        rho=cce_evolution(spin_system,generator,[],rho0,parameters.times(n),1,'final');
    else
        if segmented
            rho=rho0; scale=0;
            for k=1:numel(segments.durations{n})
                ham=segments.hamiltonians{n}{k}; dt=segments.durations{n}(k);
                pulse_prop=segments.pulses{n}{k}*expm(-1i*dt*ham);
                rho=pulse_prop*rho*pulse_prop';
                scale=scale+1+abs(dt)*norm(full(ham),2);
            end
            bound=16*spin_system.dimension^3*eps*scale;
        else
            [rho,bound]=driven_evolution(spin_system,generator,rho0,parameters.times(n),...
                                        events,starts(:,n),drives);
        end
        if controlled
            if density
                aux.abs_error_bound(n,:)=bound;
            else
                norms=cellfun(@(detector)norm(full(detector),2),detectors(:).');
                if meta.normalized, norms=norms./abs(reference); end
                aux.abs_error_bound(n,:)=bound*norms;
            end
        end
    end
    if density||~isstruct(H)
        dc=spin_system.observed_dimension; db=spin_system.bath_dimension;
        blocks=reshape(full(rho),[db dc db dc]);
        blocks=reshape(permute(blocks,[1 3 2 4]),[db*db dc*dc]);
        observed=reshape(sum(blocks(1:db+1:db*db,:),1),[dc dc]);
        if controlled, aux.density(:,:,n)=observed; end
        if density
            answer(n,:)=observed(:).';
        else
            for k=1:numel(detectors)
                if selected
                    value=exp(-1i*meta.receiver_phase)*(meta.transition.alpha'*observed*meta.transition.beta);
                else
                    detector=reshape(detectors{k},dc*db,dc*db);
                    detector=detector(1:db:end,1:db:end);
                    value=trace(detector'*observed);
                end
                if meta.normalized, value=value/reference(k); end
                answer(n,k)=value;
            end
        end
    else
        for k=1:numel(detectors)
            value=conj(detectors{k}(1))*trace(reshape(rho,spin_system.bath_dimension,[]));
            if meta.normalized, value=value/reference(k); end
            answer(n,k)=value;
        end
    end
end

end

% Compile finite events in physical-spin order without invoking a driver
function [starts,drives]=finite_controls(system,p,events)
starts=zeros(numel(events),numel(p.times));
drives=cell(1,numel(events)); local=system; local.formalism='zeeman-hilb';
for k=1:numel(events)
    e=events(k);
    if any(~ismember(e.target,system.physical_ids))
        error('cce_acquire:exterior','Finite drives require retained target spins.');
    end
    axis=e.axis;
    if ischar(axis), axis=double('xyz'==axis); end
    detuning=e.detuning;
    drive=sparse(system.dimension,system.dimension);
    for target=e.target
        for a=1:3
            names={'Lx','Ly','Lz'}; weight=e.amplitude*axis(a);
            if a==3, weight=weight+detuning; end
            if weight~=0, drive=drive+weight*cce_operator(local,names{a},target); end
        end
    end
    if strcmp(system.formalism,'zeeman-liouv')
        drive=kron(speye(system.dimension),drive)-kron(drive.',speye(system.dimension));
    end
    drives{k}=drive;
    if strcmp(e.placement,'fraction'), starts(k,:)=e.at*p.times;
    else, starts(k,:)=e.at; end
end
end

% Split every observation interval at drive starts and stops
function [rho,bound]=driven_evolution(system,L,rho,total,events,starts,drives)
stops=starts+reshape([events.duration],[],1);
scale=0; coherent=strcmp(system.formalism,'zeeman-hilb');
boundaries=unique([0;total;starts(starts>0&starts<total);stops(stops>0&stops<total)]);
if coherent, boundaries=[0;boundaries]; end
for k=1:numel(boundaries)-1
    left=boundaries(k); right=boundaries(k+1); generator=L;
    midpoint=(left+right)/2;
    active=find(starts<=midpoint&stops>midpoint&right~=left);
    for n=reshape(active,1,[]), generator=generator+drives{n}; end
    if coherent
        propagator=eye(size(generator))*expm(-1i*(right-left)*generator);
        rho=propagator*rho*propagator';
        scale=scale+1+abs(right-left)*norm(full(generator),2);
    else
        rho=cce_evolution(system,generator,[],rho,right-left,1,'final');
    end
end
bound=16*system.dimension^3*eps*scale;
end

% Propagate both conditional branches with their measured spectral residuals
function [answer,aux]=projected_fid(system,parameters,H,reference,meta)

% Diagonalize each bath Hamiltonian once in its original numerical convention
[va,ea]=eig(full((H.left+H.left')/2),'vector');
[vb,eb]=eig(full((H.right+H.right')/2),'vector');
ra=norm(H.left*va-va.*ea.','fro'); rb=norm(H.right*vb-vb.*eb.','fro');
oa=norm(va'*va-eye(size(H.left)),'fro'); ob=norm(vb'*vb-eye(size(H.right)),'fro');
residual=max(ra,rb); orthogonality=max(oa,ob);
rho=1;
for k=1:numel(system.bath_rho), rho=kron(rho,system.bath_rho{k}); end
answer=complex(zeros(numel(parameters.times),1)); bounds=zeros(size(answer));

% Keep the branch multiplication order and forward bound of the reference
for n=1:numel(parameters.times)
    time=parameters.times(n);
    ua=eye(system.bath_dimension); ub=ua;
    if time~=0
        ua=((va.*exp(-1i*time*ea).')*va')*ua;
        ub=((vb.*exp(-1i*time*eb).')*vb')*ub;
    end
    answer(n)=trace(ua*rho*ub');
    operations=12*system.bath_dimension^3; gamma=operations*eps/(1-operations*eps);
    bounds(n)=2*time*(residual+H.hamiltonian_abs_error_bound)+4*orthogonality+...
              gamma*(1+time*H.hamiltonian_norm);
end
if ~meta.normalized
    answer=answer*reference; bounds=bounds*abs(reference);
end
aux=struct('abs_error_bound',bounds,'eigen_residual',residual,...
           'orthogonality_residual',orthogonality,...
           'hamiltonian_bound',H.hamiltonian_abs_error_bound,...
           'backend',parameters.cce_aux.backend);

end

% Validate the experiment-owned time axis and retained subsystem
function grumble(spin_system,parameters,H,R,K)
if ~isstruct(spin_system)||~isfield(spin_system,'kind')||...
   ~strcmp(spin_system.kind,'cce_subsystem')||~isstruct(parameters)||...
   ~isfield(parameters,'times')||~isnumeric(parameters.times)||...
   ~isreal(parameters.times)||~isvector(parameters.times)||...
   any(~isfinite(parameters.times))||any(parameters.times<0)||any(diff(parameters.times)<0)
    error('cce_acquire:input','A native retained subsystem and increasing nonnegative times are required.');
end
if isfield(parameters,'sequence')
    error('cce_acquire:sequence','Use finite_drive, refocusing_train, or hahn_echo, not parameters.sequence.');
end
if any(isfield(parameters,{'refocusing','train'}))&&~strcmp(spin_system.formalism,'zeeman-hilb')
    error('cce_acquire:controls','Ideal refocusing requires generalized Hilbert dynamics.');
end
if ~(isnumeric(H)||isstruct(H))||~isnumeric(R)||~isnumeric(K)
    error('cce_acquire:generators','Supply conventional or tagged projected generators.');
end
end


