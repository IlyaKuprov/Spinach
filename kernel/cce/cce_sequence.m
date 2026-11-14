% Compile an explicit whole-sequence-time pulse schedule. Syntax:
%
%                   schedule=cce_sequence(model,parameters)
%
% times is a strictly increasing nonnegative seconds row starting at zero.
% sequence is an event structure array (empty means FID). An ideal event
% has type='ideal', target (physical IDs), placement='fraction'/'absolute',
% at, and either a Hermitian generator operator, applied as the rotation
% exp(-i*angle*operator) with angle defaulting to 1, or axis plus angle. axis is x/y/z
% or a real normalized three-vector. Simultaneous events retain input order.
% A selective central-register operator acts in central_ids product order;
% axis rotations are full physical-spin rotations, not qubit rotations.
% Finite controls are deliberately delegated to the cce_controls extension.
%
% Inputs:
%   model - prepared register and bath IDs.
%   parameters - times, sequence ([] for FID), method and central-state
%   declarations.
%
% Outputs:
%   schedule - sample times, compiled events and conditional/generalized
%   control metadata.
%
% See also cce_controls, cce_cluster.
%
% ilya.kuprov@weizmann.ac.il

function schedule=cce_sequence(model,parameters)

% Compile named, finite and DEER controls through the shared segment schema
advanced=isstruct(parameters.sequence)&&isfield(parameters.sequence,'name');
if isstruct(parameters.sequence)&&isfield(parameters.sequence,'type')
    advanced=advanced||any(strcmp({parameters.sequence.type},'finite'));
end
if isfield(parameters,'deer')&&isfield(parameters.deer,'enabled')
    advanced=advanced||parameters.deer.enabled;
end
if isfield(parameters,'meanfield')&&strcmp(parameters.meanfield,'static')&&...
   isstruct(parameters.sequence)&&isfield(parameters.sequence,'target')
    for event=reshape(parameters.sequence,1,[])
        if isnumeric(event.target), advanced=advanced||any(ismember(event.target,model.bath_ids)); end
    end
end
if advanced
    segments=cce_controls(model,[],[],[],parameters); schedule=segments.schedule;
    schedule.controls=strcmp(parameters.method,'gcce');
    if ~schedule.controls
        copy=parameters; copy.sequence=schedule.events;
        if any(strcmp({schedule.events.type},'finite'))
            error('cce_sequence:conditional_drive','Finite controls require gCCE.');
        end
        copy.meanfield='none';
        if isfield(copy,'deer'), copy.deer.enabled=false; end
        schedule=cce_sequence(model,copy);
    end
    return;
end

% Validate times and event syntax before producing any propagator schedule
grumble(model,parameters);

% Resolve event positions separately at every whole-sequence duration
events=parameters.sequence; nt=numel(parameters.times);
positions=zeros(numel(events),nt);
for j=1:numel(events)
    switch events(j).placement
        case 'fraction'
            positions(j,:)=events(j).at*parameters.times;
        case 'absolute'
            positions(j,:)=events(j).at;
    end
end

% Store stable ordering including the well-defined zero-duration sequence
order=cell(1,nt);
for k=1:nt
    [~,order{k}]=sort(positions(:,k),'ascend');
end
schedule=struct('times',parameters.times,'events',events,...
                'positions',positions,'order',{order},...
                'time_convention','whole_sequence_seconds',...
                'tie_rule','input_order');

end

% Validate physical targets and mutually exclusive operator definitions
function grumble(model,parameters)
if ~isstruct(model)||~isfield(model,'physical_ids')||~isstruct(parameters)||...
   ~all(isfield(parameters,{'times','sequence'}))
    error('Canonical model, times and sequence are required.');
end
t=parameters.times;
if ~isa(t,'double')||~isrow(t)||isempty(t)||any(~isfinite(t))||...
   t(1)~=0||any(diff(t)<=0)
    error('times must start at zero and increase strictly in seconds.');
end
events=parameters.sequence;
if ~isempty(events)&&~isstruct(events)
    error('sequence must be an event structure array or empty.');
end
for j=1:numel(events)
    event=events(j);
    known={'type','target','placement','at','operator','axis','angle','duration','segments'};
    if ~isempty(setdiff(fieldnames(event),known))
        error('Unknown sequence event field.');
    end
    if (isfield(event,'duration')&&~isempty(event.duration)&&event.duration~=0)||...
       (isfield(event,'segments')&&~isempty(event.segments))
        error('cce_sequence:capability','Finite duration requires cce_controls.');
    end
    if ~all(isfield(event,{'type','target','placement','at'}))||...
       ~strcmp(event.type,'ideal')
        error('cce_sequence:capability','Baseline events are ideal; finite events require cce_controls.');
    end
    if ~isnumeric(event.target)||isempty(event.target)||~isrow(event.target)||...
       numel(unique(event.target))~=numel(event.target)||...
       ~all(ismember(event.target,model.physical_ids))
        error('Pulse targets must be unique physical IDs.');
    end
    if ~ismember(event.placement,{'fraction','absolute'})||...
       ~isnumeric(event.at)||~isscalar(event.at)||~isreal(event.at)||...
       ~isfinite(event.at)||event.at<0
        error('Event placement and time must be explicit and nonnegative.');
    end
    if strcmp(event.placement,'fraction')&&event.at>1
        error('Fractional event time must lie in [0,1].');
    end
    if strcmp(event.placement,'absolute')&&any(event.at>t)
        error('Absolute event lies beyond a requested total time.');
    end
    has_operator=isfield(event,'operator')&&~isempty(event.operator);
    has_axis=isfield(event,'axis')&&~isempty(event.axis);
    if has_operator==has_axis
        error('Specify exactly one of operator or axis for each event.');
    end
    if has_operator
        d=prod(model.multiplicities(event.target)); op=event.operator;
        if ~isnumeric(op)||~isequal(size(op),[d d])||any(~isfinite(op),'all')||...
           norm(op-op','fro')>d^2*eps*max(1,norm(op,'fro'))
            error('Pulse operator must be finite Hermitian in target product order.');
        end
        if isfield(event,'angle')&&~isempty(event.angle)&&(~isnumeric(event.angle)||...
           ~isscalar(event.angle)||~isreal(event.angle)||~isfinite(event.angle))
            error('Pulse operator angle must be a finite real scalar.');
        end
    else
        if ~isfield(event,'angle')||~isnumeric(event.angle)||~isscalar(event.angle)||...
           ~isreal(event.angle)||~isfinite(event.angle)
            error('Axis rotation requires an angle in radians.');
        end
        if ischar(event.axis)
            if ~ismember(event.axis,{'x','y','z'}), error('Unknown Cartesian axis.'); end
        elseif ~isnumeric(event.axis)||~isreal(event.axis)||numel(event.axis)~=3||...
               any(~isfinite(event.axis))||abs(norm(event.axis)-1)>3*eps
            error('Numeric rotation axis must be a real unit three-vector.');
        end
    end
end
end

