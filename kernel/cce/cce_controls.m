% Compile CCE controls, selective pulses and exterior mean-field updates.
% Syntax:
%
%       segments=cce_controls(model,sites,spin_config,schedule,parameters)
%
% Empty schedule compiles parameters.sequence, including the named ideal
% sequences hahn/cpmg/xy4/xy8/udd (name,count,axis,target). Explicit finite
% events have type='finite', target, placement, at, duration (seconds),
% axis, amplitude and detuning (rad/s). They are constant rotating-frame
% drives; separate events specify piecewise envelopes. Simultaneous drives
% add, and ideal pulses preserve input order. Absolute events outside a
% requested duration are omitted (a partially intersected finite pulse is
% propagated only to that duration). Fractional events scale their start,
% never their physical width. Every sample at zero is the initial state.
%
% The returned schedule, selected_sites, durations, hamiltonians, pulses
% and exterior_expectations are indexed by time point then segment. Each
% segment propagates its Hamiltonian for duration and then applies pulse.
% Exterior ideal one-site rotations update the sampled mean field exactly;
% finite or correlation-generating exterior pulses are explicitly refused. Projected
% CCE refuses observed finite drive; general controls require gCCE.
%
% Inputs:
%   model, sites, spin_config - physical model, included bath IDs and one global
%   state assignment.
%   schedule - compiled schedule, or [] to compile parameters.
%   parameters - sequence/control fields, optional DEER and meanfield_radius
%   settings; finite events specify target, placement, at, duration in
%   seconds, axis, amplitude and detuning in rad/s.
%
% Outputs:
%   segments - schedule, selected_sites, durations, hamiltonians, pulses and
%   exterior_expectations indexed by time then segment; generalized
%   propagation also returns raw with observed densities and reducer fields.
%
% Generalized state propagation and normalization use spin_config.rho0.
%
% Internal engine records are not public experiment parameters.
% See also finite_drive, hahn_echo, cce_observation.
%
% ilya.kuprov@weizmann.ac.il

function segments=cce_controls(model,sites,spin_config,schedule,parameters)

% Validate the canonical metadata and requested controls
grumble(model,sites,spin_config,parameters);

% Compile named sequences and resolve the explicit physical event targets
if isempty(schedule)
    events=parameters.sequence;
    if isstruct(events)&&isscalar(events)&&isfield(events,'name')
        events=named_events(events);
    end
    events=select_events(model,events,parameters);
    validate_events(model,events);
    bytes=3*8*numel(events)*numel(parameters.times);
    budget=Inf;
    if isfield(parameters,'numerics')&&isfield(parameters.numerics,'max_memory_bytes')
        budget=parameters.numerics.max_memory_bytes;
    end
    if ~isfinite(bytes)||bytes>budget
        error('cce_controls:memory',...
              'Pulse positions and ordering arrays exceed numerics.max_memory_bytes.');
    end
    positions=zeros(numel(events),numel(parameters.times));
    for j=1:numel(events)
        if strcmp(events(j).placement,'fraction')
            positions(j,:)=events(j).at*parameters.times;
        else
            positions(j,:)=events(j).at;
        end
    end
    schedule=struct('events',events,'positions',positions,'times',parameters.times,...
                    'time_convention','whole_sequence_seconds');
else
    if ~isstruct(schedule)||~all(isfield(schedule,{'events','positions','times'}))||...
       ~isequal(size(schedule.positions),[numel(schedule.events) numel(schedule.times)])||...
       ~isequal(schedule.times,parameters.times)
        error('cce_controls:schedule','Compiled schedule dimensions or times do not match.');
    end
    validate_events(model,schedule.events);
end

% Publish schedule and selection even when only compilation was requested
segments=struct('schedule',schedule,'selected_sites',{arrayfun(@(e) e.target,...
                schedule.events,'UniformOutput',false)},'boundary_rule','strict_abs_less_than');
if isempty(spin_config), return; end

% Refuse approximations not represented by the driven generalized propagator
if ~ismember(parameters.method,{'cce','gcce'})||~ismember(parameters.meanfield,{'none','static'})
    error('cce_controls:method','Controls require explicit cce/gcce and none/static mean field.');
end
if isfield(parameters,'second_order')&&~isempty(parameters.second_order)&&...
   (~isstruct(parameters.second_order)||~isfield(parameters.second_order,'enabled')||parameters.second_order.enabled)
    error('cce_controls:projection','Driven controls cannot include projected second-order corrections.');
end

% Reject finite observed drives even for a zero-duration requested sample
for event=reshape(schedule.events,1,[])
    if strcmp(parameters.method,'cce')&&strcmp(event.type,'finite')&&...
       any(ismember(event.target,model.observed_ids))
        error('cce_controls:manifold_drive','Finite observed drives require gCCE.');
    end
end

% Bound retained segment matrices and exterior values before allocating operators
ids=[model.observed_ids sites]; dimension=prod(model.multiplicities(ids));
events=schedule.events; segment_count=0;
for k=1:numel(schedule.times)

    % Count distinct event boundaries without retaining any dense operator
    t=schedule.times(k); starts=schedule.positions(:,k); stops=starts;
    for j=1:numel(events)
        if strcmp(events(j).type,'finite'), stops(j)=starts(j)+events(j).duration; end
    end
    segment_count=segment_count+numel(unique([0;t;starts(starts>0&starts<t);stops(stops>0&stops<t)]));
end

% Reserve sparse complex Hamiltonians, dense pulses and exterior expectations
retained=segment_count*(48*dimension^2+24*numel(model.bath_ids)+8);
dc=prod(model.multiplicities(model.observed_ids)); components=1;
if strcmp(parameters.quantity,'density'), components=dc^2; end
if strcmp(parameters.quantity,'observable')&&iscell(parameters.coil)
    components=numel(parameters.coil);
end
retained=retained+3*numel(schedule.times)*(25*components+16*dc^2);
metadata=whos('model','sites','spin_config','schedule','parameters');
if retained+sum([metadata.bytes])+16*24*dimension^2>parameters.numerics.max_memory_bytes
    error('cce_controls:memory','Retained control segments exceed numerics.max_memory_bytes.');
end

% Build the Hamiltonian metadata only after the retained-storage check
[~,meta]=cce_hamiltonian(model,sites,struct()); ids=meta.physical_ids;
segments.durations=cell(size(schedule.times)); segments.hamiltonians=cell(size(schedule.times));
segments.pulses=cell(size(schedule.times)); segments.exterior_expectations=cell(size(schedule.times));
for t=1:numel(schedule.times)

    % Split at all ideal, drive-start and drive-stop events inside this duration
    total=schedule.times(t); starts=schedule.positions(:,t); stops=starts;
    for j=1:numel(events)
        if strcmp(events(j).type,'finite'), stops(j)=starts(j)+events(j).duration; end
    end
    boundaries=unique([0;total;starts(starts>0&starts<total);stops(stops>0&stops<total)]);
    expectations=spin_config.expectations; states=spin_config.bath_rho;
    durations=zeros(1,numel(boundaries)); hamiltonians=cell(size(durations));
    pulses=cell(size(durations)); exterior_values=cell(size(durations));
    for b=1:numel(boundaries)

        % Build the left-limit Hamiltonian for the interval ending here
        if b==1, left=0; else, left=boundaries(b-1); end
        right=boundaries(b); ham_config=struct();
        if strcmp(parameters.meanfield,'static')
            outside=~ismember(spin_config.bath_ids,sites);
            if isfield(parameters,'meanfield_radius')&&isfinite(parameters.meanfield_radius)
                distance=inf(size(spin_config.bath_ids));
                for observed_id=model.observed_ids
                    delta=model.coordinates_A(spin_config.bath_ids,:)-model.coordinates_A(observed_id,:);
                    distance=min(distance,sqrt(sum(delta.^2,2)).');
                end
                outside=outside&distance<=parameters.meanfield_radius;
            end
            ham_config.exterior_ids=spin_config.bath_ids(outside);
            ham_config.exterior_expectations=expectations(outside,:);
        end
        [H,~]=cce_hamiltonian(model,sites,ham_config);
        midpoint=(left+right)/2;
        for j=1:numel(events)
            if ~strcmp(events(j).type,'finite')||midpoint<starts(j)||midpoint>=stops(j)||right==left
                continue;
            end
            event=events(j); outside=setdiff(event.target,ids);
            if ~isempty(outside)&&strcmp(parameters.meanfield,'static')
                error('cce_controls:exterior_finite','Finite exterior pulses require dynamical exterior spins.');
            end
            if strcmp(parameters.method,'cce')&&any(ismember(event.target,model.observed_ids))
                error('cce_controls:manifold_drive','Finite observed drives require gCCE.');
            end
            axis=axis_vector(event.axis);
            for id=reshape(intersect(event.target,ids,'stable'),1,[])
                pos=ids==id; drive=zeros(dimension);
                for a=1:3, drive=drive+axis(a)*meta.operators{pos,a}; end
                H=H+event.amplitude*drive+event.detuning*meta.operators{pos,3};
            end
        end

        % Apply endpoint ideal events in stable input order and update the exterior
        pulse=eye(dimension);
        if total>0
            indices=find(starts==right);
            for j=reshape(indices,1,[])
                event=events(j);
                if ~strcmp(event.type,'ideal'), continue; end
                exterior=intersect(event.target,setdiff(spin_config.bath_ids,sites));
                if isfield(event,'operator')&&~isempty(event.operator)

                    % Omit a one-site exterior pulse only when no exterior fields are used
                    if isscalar(event.target)&&~ismember(event.target,ids)&&...
                       strcmp(parameters.meanfield,'none')
                        continue;
                    end
                    if ~isempty(exterior)
                        error('cce_controls:exterior_correlations','Exterior operator pulses require explicit dynamical sites.');
                    end
                    P=embed_operator(event.operator,event.target,ids,model.multiplicities);
                    scale=1; if isfield(event,'angle')&&~isempty(event.angle), scale=event.angle; end
                    pulse=expm(-1i*scale*P)*pulse;
                else
                    axis=axis_vector(event.axis); generator=zeros(dimension);
                    for id=reshape(intersect(event.target,ids,'stable'),1,[])
                        pos=ids==id;
                        for a=1:3, generator=generator+axis(a)*meta.operators{pos,a}; end
                    end
                    pulse=expm(-1i*event.angle*generator)*pulse;
                    for id=reshape(exterior,1,[])
                        pos=find(spin_config.bath_ids==id); spin=pauli(model.multiplicities(id)); x=spin.x; y=spin.y; z=spin.z;
                        ops={x,y,z}; P=expm(-1i*event.angle*(axis(1)*x+axis(2)*y+axis(3)*z));
                        states{pos}=P*states{pos}*P';
                        for a=1:3, expectations(pos,a)=real(trace(states{pos}*ops{a})); end
                    end
                end
            end
        end

        % Store the interval and post-pulse exterior state without hidden averaging
        durations(b)=right-left; hamiltonians{b}=H; pulses{b}=pulse;
        exterior_values{b}=expectations;
    end
    segments.durations{t}=durations; segments.hamiltonians{t}=hamiltonians;
    segments.pulses{t}=pulses; segments.exterior_expectations{t}=exterior_values;
end

% Propagate the compiled generalized Hamiltonians when requested by the engine
if strcmp(parameters.method,'gcce')
    bath=1;
    for id=reshape(sites,1,[]), bath=kron(bath,spin_config.bath_rho{spin_config.bath_ids==id}); end
    initial=kron(spin_config.rho0,bath); dc=size(spin_config.rho0,1); db=size(bath,1);
    density=zeros(dc,dc,numel(schedule.times)); bounds=zeros(numel(schedule.times),1);
    for t=1:numel(schedule.times)
        state=initial; scale=0;
        for k=1:numel(segments.durations{t})
            H=segments.hamiltonians{t}{k}; dt=segments.durations{t}(k);
            U=segments.pulses{t}{k}*expm(-1i*dt*H); state=U*state*U';
            scale=scale+1+abs(dt)*norm(full(H),2);
        end
        for a=1:dc
            for b=1:dc
                density(a,b,t)=trace(state((a-1)*db+(1:db),(b-1)*db+(1:db)));
            end
        end
        bounds(t)=16*dimension^3*eps*scale;
    end
    observe=parameters; observe.rho0=spin_config.rho0;
    if ~isfield(observe,'normalized'), observe.normalized=false; end
    value=observe_density(density,observe);
    scale=ones(1,size(value,2));
    if strcmp(parameters.quantity,'observable')
        operators=parameters.coil;
        if ~iscell(operators), operators={operators}; end
        scale=cellfun(@(operator) norm(operator,2),operators(:).');
    end
    if observe.normalized
        plain=observe; plain.normalized=false;
        scale=scale./abs(observe_density(spin_config.rho0,plain));
    end
    bounds=bounds.*scale;
    segments.raw=struct('cluster_id',uint64(0),'value',value,'valid',isfinite(value),...
                        'abs_error_bound',bounds,'density',density,'sites',sites,...
                        'state_id',spin_config.state_id,'configuration_id',spin_config.configuration_id,...
                        'eigen_residual',[],'unitarity_residual',[],...
                        'hamiltonian_backend','direct_controls');
end

end

% Generate explicit timestamps for documented named ideal sequences
function events=named_events(spec)
if ~all(isfield(spec,{'name','count','axis','target'}))||...
   ~isscalar(spec.count)||~isreal(spec.count)||~isfinite(spec.count)||...
   spec.count<1||spec.count~=fix(spec.count)
    error('cce_controls:named','Named sequences require name,count,axis,target.');
end
n=spec.count; axes=repmat({spec.axis},1,n);
switch lower(spec.name)
    case 'hahn'
        if n~=1, error('Hahn count must be one.'); end
        times=1/2;
    case 'cpmg'
        times=((1:n)-1/2)/n;
    case {'xy4','xy8'}
        if strcmpi(spec.name,'xy4'), pattern='xyxy'; else, pattern='xyxyyxyx'; end
        if mod(n,numel(pattern))~=0, error('XY count must contain complete cycles.'); end
        times=((1:n)-1/2)/n; axes=num2cell(repmat(pattern,1,n/numel(pattern)));
    case 'udd'
        times=sin(pi*(1:n)/(2*n+2)).^2;
    otherwise
        error('cce_controls:named','Unknown named pulse sequence.');
end
events=repmat(struct('type','ideal','target',spec.target,'placement','fraction',...
                    'at',0,'axis','x','angle',pi),1,n);
for k=1:n, events(k).at=times(k); events(k).axis=axes{k}; end
end

% Apply DEER selection using canonical lab-frame Zeeman frequencies in rad/s
function events=select_events(model,events,parameters)
if ~isfield(parameters,'deer')||~isfield(parameters.deer,'enabled')||~parameters.deer.enabled
    return;
end
spec=parameters.deer;
if ~all(isfield(spec,{'detuning_center','halfwidth'}))||...
   ~isscalar(spec.detuning_center)||~isreal(spec.detuning_center)||~isfinite(spec.detuning_center)||...
   ~isscalar(spec.halfwidth)||~isreal(spec.halfwidth)||~isfinite(spec.halfwidth)||spec.halfwidth<0||...
   (~isfield(spec,'sites')&&~isfield(spec,'species'))
    error('cce_controls:deer','DEER requires sites/species and a finite frequency window.');
end
candidates=model.bath_ids;
if isfield(spec,'sites')
    if any(~ismember(spec.sites,model.bath_ids))||numel(unique(spec.sites))~=numel(spec.sites)
        error('cce_controls:deer','DEER sites must be distinct bath IDs.');
    end
    candidates=intersect(candidates,spec.sites,'stable');
end
if isfield(spec,'species')
    candidates=candidates(ismember(model.spin_system.comp.isotopes(candidates),cellstr(spec.species)));
end
frequencies=zeros(size(candidates));
for k=1:numel(candidates)

    % Reject transverse canonical fields rather than silently secularizing them
    id=candidates(k); tensor=model.zeeman{id};
    if norm(tensor(3,1:2))>3*eps*norm(tensor(3,:))
        error('cce_controls:deer_axial',...
              'DEER frequency selection requires an axial Zeeman field at site %d.',id);
    end
    frequencies(k)=tensor(3,3);
end
selected=candidates(abs(frequencies-spec.detuning_center)<spec.halfwidth);
if isempty(selected)&&(~isfield(spec,'allow_empty')||~isequal(spec.allow_empty,true))
    error('cce_controls:empty_deer','DEER window selected no sites; set allow_empty explicitly.');
end
for j=1:numel(events)
    if ischar(events(j).target)&&strcmp(events(j).target,'deer')
        events(j).target=selected;
    end
end
end

% Embed an arbitrary selective observed-spin system operator in the retained tensor order
function operator=embed_operator(source,targets,ids,multiplicities)
if any(~ismember(targets,ids)), error('Selective operator targets must be retained.'); end
positions=arrayfun(@(id) find(ids==id),targets);
rest=setdiff(1:numel(ids),positions,'stable'); order=[positions rest];
dims=multiplicities(ids); d=prod(dims);
base=kron(source,eye(prod(dims(rest)))); map=zeros(d,1);
for k=0:d-1
    digits=zeros(1,numel(dims)); remainder=k;
    for j=numel(dims):-1:1, digits(j)=mod(remainder,dims(j)); remainder=floor(remainder/dims(j)); end
    index=0;
    for j=order, index=index*dims(j)+digits(j); end
    map(k+1)=index+1;
end
operator=base(map,map);
end

% Convert validated axes into real Cartesian vectors
function axis=axis_vector(value)
if ischar(value)
    axis=zeros(1,3); axis('xyz'==value)=1;
else
    axis=reshape(value,1,3);
end
end

% Validate events after named generation and frequency-window mapping
function validate_events(model,events)
if isempty(events), return; end
if ~isstruct(events), error('Control events must be structures.'); end
for event=reshape(events,1,[])
    if ~all(isfield(event,{'type','target','placement','at'}))||...
       ~ismember(event.type,{'ideal','finite'})||~ismember(event.placement,{'fraction','absolute'})||...
       ~isnumeric(event.target)||~isreal(event.target)||numel(unique(event.target))~=numel(event.target)||...
       any(~ismember(event.target,[model.observed_ids model.bath_ids]))||...
       ~isscalar(event.at)||~isreal(event.at)||~isfinite(event.at)||event.at<0||...
       (strcmp(event.placement,'fraction')&&event.at>1)
        error('cce_controls:event','Invalid type, target or placement.');
    end
    if strcmp(event.type,'finite')
        if ~all(isfield(event,{'duration','axis','amplitude','detuning'}))||...
           ~isscalar(event.duration)||~isreal(event.duration)||~isfinite(event.duration)||event.duration<=0||...
           ~isscalar(event.amplitude)||~isreal(event.amplitude)||~isfinite(event.amplitude)||...
           ~isscalar(event.detuning)||~isreal(event.detuning)||~isfinite(event.detuning)
            error('cce_controls:finite','Finite pulses require duration,axis,amplitude,detuning.');
        end
    else

        % Operator events default to angle one; axis events require an angle
        has_operator=isfield(event,'operator')&&~isempty(event.operator);
        has_angle=isfield(event,'angle')&&~isempty(event.angle);
        if (~has_operator&&~has_angle)||(has_angle&&(~isnumeric(event.angle)||...
           ~isscalar(event.angle)||~isreal(event.angle)||~isfinite(event.angle)))
            error('cce_controls:angle','Ideal pulse angle must be finite and real.');
        end
    end
    if strcmp(event.type,'ideal')&&isfield(event,'operator')&&~isempty(event.operator)
        d=prod(model.multiplicities(event.target)); op=event.operator;
        if ~isnumeric(op)||~isequal(size(op),[d d])||any(~isfinite(op),'all')||...
           norm(op-op','fro')>d^2*eps*norm(op,'fro')
            error('cce_controls:operator','Pulse operator must be finite and Hermitian.');
        end
    elseif ~isfield(event,'axis')||...
           ~(ischar(event.axis)&&ismember(event.axis,{'x','y','z'}))&&...
           ~(isnumeric(event.axis)&&isreal(event.axis)&&numel(event.axis)==3&&...
             all(isfinite(event.axis))&&abs(norm(event.axis)-1)<=3*eps)
        error('cce_controls:axis','Pulse axis must be x/y/z or a real unit vector.');
    end
end
end

% Extract density components and retain the explicit reference normalisation
function observed=observe_density(rho,parameters)

% Validate the density stack and requested observable
detection_grumble(rho,parameters);

% Extract each time point in an explicitly documented element order
d=size(rho,1); nt=size(rho,3);
switch parameters.quantity
    case 'density'
        observed=reshape(rho,d*d,nt).';
    case 'signal'
        observed=complex(zeros(nt,1));
        for k=1:nt
            observed(k)=parameters.transition.alpha'*rho(:,:,k)*parameters.transition.beta;
        end
    case 'observable'
        coils=parameters.coil;
        if ~iscell(coils), coils={coils}; end
        observed=complex(zeros(nt,numel(coils)));
        for k=1:nt
            for j=1:numel(coils)
                observed(k,j)=trace(coils{j}'*rho(:,:,k));
            end
        end
end

% Normalize only against the explicitly supplied initial observed state
if isfield(parameters,'normalized')&&parameters.normalized
    if strcmp(parameters.quantity,'density')
        error('cce_observe:normalization','A density matrix is not normalized elementwise.');
    end
    reference_parameters=parameters; reference_parameters.normalized=false;
    reference=observe_density(parameters.rho0,reference_parameters);
    bound=size(rho,1)^2*eps*max(1,norm(parameters.rho0,'fro'));

    % Propagate density uncertainty in the units of each observable
    if strcmp(parameters.quantity,'observable')
        bound=bound*cellfun(@(op) norm(op,'fro'),coils(:).');
    end
    if any(abs(reference)<=bound)
        error('cce_observe:initial_zero','Initial observable is unresolved or zero.');
    end
    observed=observed./reference;
end

end

% Validate dimensions, finite values and observation definitions
function detection_grumble(rho,parameters)
if ~isnumeric(rho)||ndims(rho)>3||size(rho,1)~=size(rho,2)||...
   isempty(rho)||any(~isfinite(rho),'all')
    error('rho must be a finite square matrix or stack.');
end
if ~isstruct(parameters)||~isscalar(parameters)||~isfield(parameters,'quantity')||...
   ~ismember(parameters.quantity,{'density','signal','observable'})
    error('quantity must be density, signal, or observable.');
end
if strcmp(parameters.quantity,'signal')
    if ~isfield(parameters,'transition')||...
       ~all(isfield(parameters.transition,{'alpha','beta'}))
        error('Signal detection requires parameters.transition.alpha and parameters.transition.beta.');
    end
    for name={'alpha','beta'}
        state=parameters.transition.(name{1});
        if ~isnumeric(state)||~isequal(size(state),[size(rho,1) 1])||...
           any(~isfinite(state))||abs(norm(state)-1)>size(rho,1)*eps
            error('Observed-spin observation states must be normalized columns.');
        end
    end
    if abs(parameters.transition.alpha'*parameters.transition.beta)>size(rho,1)*eps
        error('Observed-spin observation states must be orthogonal.');
    end
end
if strcmp(parameters.quantity,'observable')
    if ~isfield(parameters,'coil'), error('coil must be supplied.'); end
    coils=parameters.coil;
    if ~iscell(coils), coils={coils}; end
    if isempty(coils), error('coils must be nonempty.'); end
    for k=1:numel(coils)
        if ~isnumeric(coils{k})||~isequal(size(coils{k}),size(rho,[1 2]))||...
           any(~isfinite(coils{k}),'all')
            error('Each receiver coil must match the observed-spin dimension.');
        end
    end
end
if isfield(parameters,'normalized')
    if ~islogical(parameters.normalized)||~isscalar(parameters.normalized)
        error('normalized must be a logical scalar.');
    end
    if parameters.normalized&&~isfield(parameters,'rho0')
        error('Normalization requires rho0.');
    end
end
end

% Check canonical metadata and physical site IDs
function grumble(model,sites,spin_config,parameters)
if ~isstruct(model)||~isscalar(model)||~all(isfield(model,{'observed_ids','bath_ids','multiplicities'}))||...
   ~isnumeric(sites)||~isreal(sites)||any(~ismember(sites,model.bath_ids))||numel(unique(sites))~=numel(sites)||...
   ~isstruct(parameters)||~all(isfield(parameters,{'times','sequence','method','meanfield'}))||...
   ~isnumeric(parameters.times)||~isreal(parameters.times)||~isrow(parameters.times)||...
   isempty(parameters.times)||parameters.times(1)~=0||any(~isfinite(parameters.times))||...
   any(diff(parameters.times)<=0)
    error('cce_controls:input','Canonical model, sites and increasing times starting at zero required.');
end
if ~isempty(spin_config)&&(~isstruct(spin_config)||~all(isfield(spin_config,{'bath_ids','bath_rho','expectations'})))
    error('cce_controls:spin_config','A canonical spin-configuration structure is required.');
end
if isfield(parameters,'initial_rho')
    error('cce:spin_config','Use spin_config.rho0 for the observed density.');
end
if ~isempty(spin_config)
    known={'rho0','bath_ids','bath_rho','expectations','state_id',...
           'configuration_id','weight','assignment','id'};
    if ~isscalar(spin_config)||~isempty(setdiff(fieldnames(spin_config),known))
        error('cce:spin_config','Unsupported spin-configuration field.');
    end
end
if isfield(parameters,'meanfield_radius')
    radius=parameters.meanfield_radius;
    if ~isnumeric(radius)||~isreal(radius)||~isscalar(radius)||isnan(radius)||radius<0
        error('cce:parameters','meanfield_radius must be a nonnegative real scalar.');
    end
    if ~isempty(spin_config)&&isfinite(radius)
        exterior=setdiff(spin_config.bath_ids,sites);
        if ~isempty(exterior)&&(~isfield(model,'coordinates_A')||...
           ~isnumeric(model.coordinates_A)||~isreal(model.coordinates_A)||...
           size(model.coordinates_A,2)~=3||...
           size(model.coordinates_A,1)<max([model.observed_ids exterior])||...
           any(~isfinite(model.coordinates_A([model.observed_ids exterior],:)),'all'))
            error('cce:coordinates','Finite observed and exterior coordinates are required.');
        end
    end
end
end


