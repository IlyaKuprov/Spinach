% Compile experiment-owned ideal and finite controls in a native subsystem.
% Syntax:
%
%       segments=cce_segments(system,parameters,drift,events)
%
% Internal service: events are compiled by an experiment, not a public input
% interface. Each sample retains durations, angular-frequency Hamiltonians,
% endpoint pulse matrices, and post-pulse exterior spin expectations. Drift
% remains active during finite pulses. Zero elapsed time retains the initial
% state; absolute pulses outside an endpoint are omitted. Static exterior
% Cartesian rotations update subsequent fields; exterior operator pulses
% and finite exterior drives are refused. Memory is checked before matrices.
%
% Parameters:
%
%   system     - retained subsystem returned by cce_generators
%   parameters - experiment parameters with a seconds-valued times row
%                and cce controls
%   drift      - retained drift Hamiltonian in rad/s
%   events     - internal finite/ideal event structure row
%
% Outputs:
%
%   segments - timing schedule and one cell per sample for interval durations
%              in seconds, Hamiltonians in rad/s, endpoint pulses, and
%              exterior spin expectations
%
% ilya.kuprov@weizmann.ac.il

function segments=cce_segments(system,parameters,drift,events)

% Validate internal events and construct their physical timing table
grumble(system,parameters,drift,events);
positions=zeros(numel(events),numel(parameters.times));
for k=1:numel(events)
    if strcmp(events(k).placement,'fraction'), positions(k,:)=events(k).at*parameters.times;
    else
        positions(k,:)=events(k).at;
    end
end
schedule=struct('events',events,'positions',positions,'times',parameters.times);
segments=struct('schedule',schedule);

% Bound retained segment matrices and exterior values before allocating operators
sites=system.cce.sites;
spin_config=system.cce.configuration;
ids=system.cce.ids;
dimension=system.cce.dimension;
events=schedule.events;
segment_count=0;
for k=1:numel(schedule.times)

    % Count distinct event boundaries without retaining any dense operator
    t=schedule.times(k);
    starts=schedule.positions(:,k);
    stops=starts;
    for j=1:numel(events)
        if strcmp(events(j).type,'finite')
            stops(j)=starts(j)+events(j).duration;
        end
    end
    segment_count=segment_count+numel(unique([0;t;starts(starts>0&starts<t);stops(stops>0&stops<t)]));

end

% Reserve sparse complex Hamiltonians, dense pulses, and exterior expectations
retained=segment_count*(48*dimension^2+24*numel(parameters.cce.bath_spins)+8);
dc=prod(system.comp.mults(parameters.observed));
components=1;
if isfield(parameters,'output')&&strcmp(parameters.output,'density')
    components=dc^2;
end
if isfield(parameters,'output')&&strcmp(parameters.output,'observable')&&...
   isfield(parameters,'coil')&&iscell(parameters.coil)
    components=numel(parameters.coil);
end
retained=retained+3*numel(schedule.times)*(25*components+16*dc^2);
metadata=whos('system','schedule','parameters');
if retained+sum([metadata.bytes])+16*24*dimension^2>system.tols.cce_mem
    error('cce_segments:memory','Retained control segments exceed sys.tols.cce_mem.');
end

% Build the Hamiltonian metadata only after the retained-storage check
operators=system.cce.operators;
segments.durations=cell(size(schedule.times));
segments.hamiltonians=cell(size(schedule.times));
segments.pulses=cell(size(schedule.times));
segments.exterior_expectations=cell(size(schedule.times));
for t=1:numel(schedule.times)

    % Split at all ideal, drive-start, and drive-stop events inside this duration
    total=schedule.times(t);
    starts=schedule.positions(:,t);
    stops=starts;
    for j=1:numel(events)
        if strcmp(events(j).type,'finite')
            stops(j)=starts(j)+events(j).duration;
        end
    end
    boundaries=unique([0;total;starts(starts>0&starts<total);stops(stops>0&stops<total)]);
    expectations=spin_config.expectations;
    states=spin_config.bath_rho;
    durations=zeros(1,numel(boundaries));
    hamiltonians=cell(size(durations));
    pulses=cell(size(durations));
    exterior_values=cell(size(durations));
    for b=1:numel(boundaries)

        % Build the left-limit Hamiltonian for the interval ending here
        if b==1
            left=0;
        else
            left=boundaries(b-1);
        end
        right=boundaries(b);
        ham_config=spin_config;
        if strcmp(parameters.cce.meanfield,'static')
            ham_config.bath_rho=states;
            for field={'expectations','exterior_ids','exterior_expectations'}
                if isfield(ham_config,field{1})
                    ham_config=rmfield(ham_config,field{1});
                end
            end
            [current,~]=cce_subsystem(system,parameters,sites,ham_config,[]);
            [H,~,~]=cce_generators(current,parameters,[]);
        else
            H=drift;
        end
        midpoint=(left+right)/2;
        for j=1:numel(events)
            if ~strcmp(events(j).type,'finite')||midpoint<starts(j)||midpoint>=stops(j)||right==left
                continue;
            end
            event=events(j);
            axis=axis_vector(event.axis);
            for id=reshape(intersect(event.target,ids,'stable'),1,[])
                pos=ids==id;
                drive=zeros(dimension);
                for a=1:3
                    drive=drive+axis(a)*operators{ids==id,a};
                end
                H=H+event.amplitude*drive+event.detuning*operators{pos,3};
            end
        end

        % Apply endpoint ideal events in stable input order and update the exterior
        pulse=eye(dimension);
        if total>0
            indices=find(starts==right);
            for j=reshape(indices,1,[])
                event=events(j);
                if ~strcmp(event.type,'ideal')
                    continue;
                end
                exterior=intersect(event.target,setdiff(spin_config.bath_ids,sites));
                if isfield(event,'operator')&&~isempty(event.operator)

                    % Omit a one-site exterior pulse only when no exterior fields are used
                    if isscalar(event.target)&&~ismember(event.target,ids)&&...
                       strcmp(parameters.cce.meanfield,'none')
                        continue;
                    end

                    % Embed the selective operator in the retained tensor order
                    slots=arrayfun(@(id) find(ids==id),event.target);
                    rest=setdiff(1:numel(ids),slots,'stable');
                    order=[slots rest];
                    dims=system.comp.mults(ids);
                    d=prod(dims);
                    base=kron(event.operator,eye(prod(dims(rest))));
                    map=zeros(d,1);
                    for k=0:d-1
                        digits=zeros(1,numel(dims));
                        remainder=k;
                        for q=numel(dims):-1:1
                            digits(q)=mod(remainder,dims(q));
                            remainder=floor(remainder/dims(q));
                        end
                        index=0;
                        for q=order
                            index=index*dims(q)+digits(q);
                        end
                        map(k+1)=index+1;
                    end
                    P=base(map,map);
                    scale=1;
                    if isfield(event,'angle')&&~isempty(event.angle)
                        scale=event.angle;
                    end
                    pulse=expm(-1i*scale*P)*pulse;
                else
                    axis=axis_vector(event.axis);
                    generator=zeros(dimension);
                    for id=reshape(intersect(event.target,ids,'stable'),1,[])
                        for a=1:3
                            generator=generator+axis(a)*operators{ids==id,a};
                        end
                    end
                    pulse=expm(-1i*event.angle*generator)*pulse;
                    for id=reshape(exterior,1,[])
                        pos=find(spin_config.bath_ids==id);
                        spin=pauli(system.comp.mults(id));
                        x=spin.x;
                        y=spin.y;
                        z=spin.z;
                        ops={x,y,z};
                        P=expm(-1i*event.angle*(axis(1)*x+axis(2)*y+axis(3)*z));
                        states{pos}=P*states{pos}*P';
                        for a=1:3
                            expectations(pos,a)=real(trace(states{pos}*ops{a}));
                        end
                    end

                end
            end
        end

        % Store the interval and post-pulse exterior state without hidden averaging
        durations(b)=right-left;
        hamiltonians{b}=H;
        pulses{b}=pulse;
        exterior_values{b}=expectations;

    end
    segments.durations{t}=durations;
    segments.hamiltonians{t}=hamiltonians;
    segments.pulses{t}=pulses;
    segments.exterior_expectations{t}=exterior_values;

end

end

% Convert validated axes into real Cartesian vectors
function axis=axis_vector(value)
if ischar(value)
    axis=zeros(1,3);
    axis('xyz'==value)=1;
else
    axis=reshape(value,1,3);
end
end

% Validate experiment-compiled events before allocating retained matrices
function grumble(system,parameters,drift,events)
if isfield(parameters,'sequence')
    error('cce_segments:sequence','Use finite_drive or refocusing_train native controls, not parameters.sequence.');
end
if ~isstruct(system)||~isfield(system,'cce')||~isfield(system.cce,'sites')||...
   ~isstruct(parameters)||~isfield(parameters,'times')||~isnumeric(drift)
    error('cce_segments:input','Supply a native subsystem, parameters, and numeric drift.');
end
if ~isstruct(events)||~isrow(events)
    error('cce_segments:events','The experiment must supply an event structure row.');
end
for event=events
    if ~all(isfield(event,{'type','target','placement','at'}))||...
       ~ismember(event.type,{'ideal','finite'})||...
       ~isnumeric(event.target)||~isrow(event.target)||isempty(event.target)||...
       numel(unique(event.target))~=numel(event.target)||...
       any(~ismember(event.target,1:system.comp.nspins))||...
       ~ismember(event.placement,{'absolute','fraction'})||...
       ~isnumeric(event.at)||~isreal(event.at)||~isscalar(event.at)||...
       ~isfinite(event.at)||event.at<0||...
       (strcmp(event.placement,'fraction')&&event.at>1)
        error('cce_segments:events','Invalid physical target or pulse timing.');
    end
    if strcmp(event.type,'finite')
        if ~all(isfield(event,{'duration','amplitude','detuning','axis'}))||...
           ~isnumeric([event.duration event.amplitude event.detuning])||...
           ~isreal([event.duration event.amplitude event.detuning])||...
           numel([event.duration event.amplitude event.detuning])~=3||...
           any(~isfinite([event.duration event.amplitude event.detuning]))||event.duration<=0
            error('cce_segments:finite','Finite pulses require positive duration and finite rad/s amplitudes.');
        end
        if strcmp(parameters.cce.meanfield,'static')&&...
           any(~ismember(event.target,system.cce.ids))
            error('cce_segments:exterior_finite','Finite exterior pulses require dynamical exterior spins.');
        end
        if strcmp(parameters.cce.approximation,'projected')&&...
           any(ismember(event.target,parameters.observed))
            error('cce_segments:manifold_drive','Finite observed drives require gCCE.');
        end
    elseif isfield(event,'operator')&&~isempty(event.operator)
        d=prod(system.comp.mults(event.target));
        if ~isnumeric(event.operator)||~isequal(size(event.operator),[d d])||...
           any(~isfinite(event.operator),'all')||...
           norm(event.operator-event.operator','fro')>d*eps*norm(event.operator,'fro')
            error('cce_segments:operator','Pulse generators must be finite and Hermitian.');
        end
        outside=~ismember(event.target,system.cce.ids);
        if any(outside)&&~(isscalar(event.target)&&strcmp(parameters.cce.meanfield,'none'))
            if any(ismember(event.target(outside),system.cce.configuration.bath_ids))
                error('cce_segments:exterior_correlations','Exterior operator pulses require explicit dynamical sites.');
            end
            error('Selective operator targets must be retained.');
        end
    end
    if ~isfield(event,'operator')||isempty(event.operator)
        axis=event.axis;
        if ischar(axis)&&ismember(axis,{'x','y','z'})
            axis=double('xyz'==axis);
        end
        if ~isnumeric(axis)||~isreal(axis)||~isequal(size(axis),[1 3])||...
           any(~isfinite(axis))||abs(norm(axis)-1)>3*eps
            error('cce_segments:axis','Supply a Cartesian label or a real unit axis.');
        end
    end
    if strcmp(event.type,'ideal')&&(~isfield(event,'angle')||...
       ~isnumeric(event.angle)||~isreal(event.angle)||~isscalar(event.angle)||~isfinite(event.angle))
        error('cce_segments:angle','An ideal event requires a finite real angle.');
    end
end
end


