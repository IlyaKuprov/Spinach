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
coherent=isempty(events)&&~isstruct(H)&&strcmp(spin_system.formalism,'zeeman-hilb');
liouville=density&&isempty(events)&&~isstruct(H)&&...
          strcmp(spin_system.formalism,'zeeman-liouv');

% Preserve the stock density action and its arithmetic estimate
if liouville
    engine=spin_system.model.spin_system; engine.bas.formalism='sphten-liouv';
    engine.tols.prop_chop=0;
    scale=ones(1,size(answer,2));
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
            rho=step(engine,generator,rho0,parameters.times(n));
        end
    elseif isempty(events)
        rho=cce_evolution(spin_system,generator,[],rho0,parameters.times(n),1,'final');
    else
        rho=driven_evolution(spin_system,generator,rho0,parameters.times(n),...
                             events,starts(:,n),drives);
    end
    if density
        dc=spin_system.observed_dimension; db=spin_system.bath_dimension;
        rho=reshape(rho,dc*db,dc*db); observed=zeros(dc);
        for a=1:dc
            for b=1:dc
                observed(a,b)=trace(rho((a-1)*db+(1:db),(b-1)*db+(1:db)));
            end
        end
        answer(n,:)=observed(:).';
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
function rho=driven_evolution(system,L,rho,total,events,starts,drives)
stops=starts+reshape([events.duration],[],1);
boundaries=unique([0;total;starts(starts>0&starts<total);stops(stops>0&stops<total)]);
for k=1:numel(boundaries)-1
    left=boundaries(k); right=boundaries(k+1); generator=L;
    active=find(starts<right&stops>left);
    for n=reshape(active,1,[]), generator=generator+drives{n}; end
    rho=cce_evolution(system,generator,[],rho,right-left,1,'final');
end
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


