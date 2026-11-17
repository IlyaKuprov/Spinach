% Free evolution of a retained cluster with experiment-owned detection.
% Syntax:
%
%     [answer,aux]=free_evolution(spin_system,parameters,H,R,K)
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
% sequence is an optional finite-event array: target physical IDs, axis
% (x/y/z or a unit vector), amplitude and detuning in rad/s, duration in
% seconds, placement (absolute/fraction) and at (seconds/fraction). Drift
% and relaxation remain active during each drive. Overlapping drives add.
% Drives on exterior bath spins and projected dynamics are not supported.
% H,R,K are conventional angular-frequency generators or explicitly tagged
% projected branches from cce_generators. No pulse or receiver is implicit.
%
% GKSL relaxation uses L=H+1i*R+1i*K. Projected propagation carries a
% coherence block, not a physical density; density requires generalized CCE.
%
% ilya.kuprov@weizmann.ac.il

function [answer,aux]=free_evolution(spin_system,parameters,H,R,K)

% Preparation and detection belong to this experiment, not to its context
grumble(spin_system,parameters,H,R,K);
if ~isfield(parameters,'normalized'), parameters.normalized=true; end
if isfield(parameters,'output')&&strcmp(parameters.output,'density')
    parameters.normalized=false;
    if ~isfield(parameters,'coil'), parameters.coil=eye(spin_system.observed_dimension); end
end
[rho0,coil,reference,meta]=cce_observation(spin_system,parameters);
selected=~isfield(parameters,'coil')&&~isfield(parameters,'detection')&&isfield(meta,'transition');

% Preserve the projected spectral action for a selected-transition FID
if isstruct(H)&&isempty(H.dissipator)&&...
   ~(isfield(parameters,'output')&&strcmp(parameters.output,'density'))&&(~isfield(parameters,'sequence')||isempty(parameters.sequence))&&...
   ~isfield(parameters,'rho0')&&~isfield(parameters,'coil')&&~isfield(parameters,'detection')&&...
   ~isfield(parameters,'preparation')&&~isfield(parameters,'bath_operator')
    [answer,aux]=projected_fid(spin_system,parameters,H,reference,meta); return
end
if isstruct(H), generator=H; else, generator=H+1i*R+1i*K; end
density=isfield(parameters,'output')&&strcmp(parameters.output,'density');
if density
    if isstruct(H), error('free_evolution:density','Projected branches do not return a density matrix.'); end
    answer=complex(zeros(numel(parameters.times),spin_system.observed_dimension^2));
else
    answer=complex(zeros(numel(parameters.times),numel(reference)));
end

% Compile experiment-owned finite drives once for this retained subsystem
[events,starts,drives]=finite_controls(spin_system,parameters,H);

% Reuse the coherent spectrum and retain its measured numerical diagnostics
aux=struct();
controlled=~isempty(events)&&~isstruct(H)&&strcmp(spin_system.formalism,'zeeman-hilb');
if controlled, aux.abs_error_bound=zeros(size(answer)); end
coherent=isempty(events)&&~isstruct(H)&&strcmp(spin_system.formalism,'zeeman-hilb');
liouville=isempty(events)&&~isstruct(H)&&...
          strcmp(spin_system.formalism,'zeeman-liouv');

% Preserve the stock density action and its arithmetic estimate
if liouville
    engine=spin_system.model.spin_system; engine.bas.formalism='sphten-liouv';
    engine.tols.prop_chop=0;
    scale=ones(1,size(answer,2));
    if ~density
        if ~selected
            if ~iscell(coil), detectors={coil}; else, detectors=coil; end
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
        if ~iscell(coil), detectors={coil}; else, detectors=coil; end
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
            if isdiag(action), action=sparse(action); end
            rho=step(engine,action,rho0,parameters.times(n));
        end
    elseif isempty(events)
        rho=cce_evolution(spin_system,generator,[],rho0,parameters.times(n),1,'final');
    else
        [rho,bound]=driven_evolution(spin_system,generator,rho0,parameters.times(n),...
                             events,starts(:,n),drives);
        if controlled
            if density
                aux.abs_error_bound(n,:)=bound;
            else
                norms=cellfun(@(detector)norm(full(detector),2),detectors);
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
        if density
            answer(n,:)=observed(:).';
        else
            if ~iscell(coil), detectors={coil}; else, detectors=coil; end
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
        if ~iscell(coil), detectors={coil}; else, detectors=coil; end
        for k=1:numel(detectors)
            if isvector(rho), value=detectors{k}'*rho;
            else, value=trace(detectors{k}'*rho); end
            if meta.normalized, value=value/reference(k); end
            answer(n,k)=value;
        end
    end
end

end

% Compile finite events in physical-spin order without invoking a driver
function [events,starts,drives]=finite_controls(system,p,H)
events=struct([]); starts=zeros(0,numel(p.times)); drives={};
if ~isfield(p,'sequence')||isempty(p.sequence), return; end
if ~isstruct(p.sequence)||any(~isfield(p.sequence,{'type','target','axis',...
        'amplitude','placement','at','duration'}))
    error('free_evolution:sequence','Supply finite events with physical targets and timings.');
end
if isstruct(H)
    error('free_evolution:projection','Finite drives require generalized evolution.');
end
events=p.sequence; starts=zeros(numel(events),numel(p.times));
drives=cell(1,numel(events)); local=system; local.formalism='zeeman-hilb';
for k=1:numel(events)
    e=events(k);
    if ~strcmp(e.type,'finite')||~isnumeric(e.target)||~isrow(e.target)||...
       isempty(e.target)||numel(unique(e.target))~=numel(e.target)||...
       any(~ismember(e.target,system.model.physical_ids))
        error('free_evolution:sequence','Only finite drives on distinct physical spin IDs are accepted.');
    end
    if any(~ismember(e.target,system.physical_ids))
        error('free_evolution:exterior','Finite drives require retained target spins.');
    end
    for name={'amplitude','at','duration'}
        v=e.(name{1});
        if ~isnumeric(v)||~isreal(v)||~isscalar(v)||~isfinite(v)
            error('free_evolution:sequence','Finite drive amplitudes and timings must be real finite scalars.');
        end
    end
    if e.duration<=0||e.at<0||~ismember(e.placement,{'absolute','fraction'})||...
       (strcmp(e.placement,'fraction')&&e.at>1)
        error('free_evolution:sequence','Invalid drive duration or placement.');
    end
    axis=e.axis;
    if ischar(axis)&&ismember(axis,{'x','y','z'})
        axis=double('xyz'==axis);
    end
    if ~isnumeric(axis)||~isreal(axis)||numel(axis)~=3||...
       any(~isfinite(axis))||abs(norm(axis)-1)>3*eps
        error('free_evolution:sequence','A drive axis is x, y, z or a unit three-vector.');
    end
    detuning=0;
    if isfield(e,'detuning')&&~isempty(e.detuning), detuning=e.detuning; end
    if ~isnumeric(detuning)||~isreal(detuning)||~isscalar(detuning)||~isfinite(detuning)
        error('free_evolution:sequence','Drive detuning must be a real finite rad/s scalar.');
    end
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
    active=find(starts<right&stops>left);
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
    error('free_evolution:input','A native retained subsystem and increasing nonnegative times are required.');
end
if ~(isnumeric(H)||isstruct(H))||~isnumeric(R)||~isnumeric(K)
    error('free_evolution:generators','Supply conventional or tagged projected generators.');
end
end


