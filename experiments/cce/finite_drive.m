% Finite rectangular drive followed by free evolution in a retained cluster.
% Syntax:
%
%     [answer,aux]=finite_drive(spin_system,parameters,H,R,K)
%
% parameters.times gives elapsed times in seconds, each zero or at least
% drive.duration. parameters.drive is a scalar structure with fields:
%   spins     - unique physical observed-spin identifiers
%   axis      - 'x', 'y', 'z', or a real unit three-vector
%   amplitude - real finite drive amplitude in rad/s
%   detuning  - real finite longitudinal offset in rad/s
%   duration  - positive finite rectangular-pulse duration in seconds
% The drive starts at time zero; the complete drift stays active. Preparation,
% detection, output and normalization follow free_evolution. Projected branches
% are not supported. Pulse events are compiled inside this experiment.
%
% ilya.kuprov@weizmann.ac.il

function [answer,aux]=finite_drive(spin_system,parameters,H,R,K)

% Validate the declarative drive before compiling its finite event
grumble(spin_system,parameters,H,R,K);
drive=parameters.drive;
parameters.sequence=struct('type','finite','target',drive.spins,...
    'axis',drive.axis,'amplitude',drive.amplitude,'detuning',drive.detuning,...
    'placement','fraction','at',0,'duration',drive.duration);

% Reuse the experiment-level finite-control propagation contract
[answer,aux]=free_evolution(spin_system,parameters,H,R,K);

end

% Validate every drive field and reject competing pulse specifications
function grumble(spin_system,parameters,H,R,K)
if ~isstruct(spin_system)||~isfield(spin_system,'model')||...
   ~isstruct(parameters)||~isscalar(parameters)||~isfield(parameters,'drive')
    error('finite_drive:input','Supply a retained subsystem and drive parameters.');
end
if isfield(parameters,'sequence')&&~isempty(parameters.sequence)
    error('finite_drive:sequence','Pulse events belong to the finite_drive experiment.');
end
drive=parameters.drive;
fields={'spins','axis','amplitude','detuning','duration'};
if ~isstruct(drive)||~isscalar(drive)||any(~isfield(drive,fields))||...
   ~isempty(setdiff(fieldnames(drive),fields))
    error('finite_drive:drive','Supply exactly spins, axis, amplitude, detuning and duration.');
end
if ~isnumeric(drive.spins)||~isreal(drive.spins)||~isrow(drive.spins)||...
   isempty(drive.spins)||any(~isfinite(drive.spins))||...
   any(mod(drive.spins,1)~=0)||numel(unique(drive.spins))~=numel(drive.spins)||...
   any(~ismember(drive.spins,spin_system.model.observed_ids))
    error('finite_drive:spins','Drive spins must be unique physical observed-spin identifiers.');
end
axis=drive.axis;
if ~(ischar(axis)&&ismember(axis,{'x','y','z'}))&&...
   ~(isnumeric(axis)&&isreal(axis)&&isvector(axis)&&numel(axis)==3&&...
     all(isfinite(axis))&&abs(norm(axis)-1)<=3*eps)
    error('finite_drive:axis','Drive axis must be x, y, z or a real unit three-vector.');
end
for name={'amplitude','detuning','duration'}
    value=drive.(name{1});
    if ~isnumeric(value)||~isreal(value)||~isscalar(value)||~isfinite(value)
        error('finite_drive:scalar','Drive amplitudes, offsets and durations must be real finite scalars.');
    end
end
if drive.duration<=0
    error('finite_drive:duration','Drive duration must be positive.');
end
if ~isfield(parameters,'times')||~isnumeric(parameters.times)||...
   ~isreal(parameters.times)||~isvector(parameters.times)||isempty(parameters.times)||...
   any(~isfinite(parameters.times))||any(parameters.times<0)||...
   any(diff(parameters.times)<0)||...
   any(parameters.times>0&parameters.times<drive.duration)
    error('finite_drive:times','Times must increase, and be zero or at least the drive duration.');
end
if ~isnumeric(H)||~isnumeric(R)||~isnumeric(K)
    error('finite_drive:generators','Finite drives require conventional numeric generators.');
end
end


