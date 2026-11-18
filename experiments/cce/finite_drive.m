% Finite rectangular or piecewise drive in a retained cluster.
% Syntax:
%
%     [answer,aux]=finite_drive(spin_system,parameters,H,R,K)
%
% parameters.times gives increasing elapsed seconds, including during a drive.
% parameters.drive is a nonempty structure row, each element with fields:
%   spins     - unique physical observed-spin identifiers
%   axis      - 'x', 'y', 'z', or a real unit three-vector
%   amplitude - real finite amplitude row in rad/s
%   detuning  - real finite longitudinal offset in rad/s
%   duration  - matching positive duration row in seconds
% Optional time specifies the start in seconds; otherwise the start is zero.
% Alternatively fraction specifies one start as a fraction of each elapsed
% time, for a single finite pulse whose duration remains in seconds.
% Equal-length amplitude/duration rows specify consecutive waveform slices.
% A matching time row gives explicit starts, including overlapping slices.
% Optional refocusing is a cell row of ideal-pulse descriptors (spins,
% operator or axis/angle, and time or fraction). The complete drift stays active. Preparation,
% detection, output and normalization follow free_evolution. Projected branches
% are not supported. Pulse events are compiled inside this experiment.
% Mixed ideal/finite controls require generalized Hilbert dynamics.
%
% ilya.kuprov@weizmann.ac.il

function [answer,aux]=finite_drive(spin_system,parameters,H,R,K)

% Validate the declarative drive before compiling its finite event
grumble(spin_system,parameters,H,R,K);

% Reuse the experiment-level finite-control propagation contract
[answer,aux]=cce_acquire(spin_system,parameters,H,R,K);

end

% Validate every drive field and reject competing pulse specifications
function grumble(spin_system,parameters,H,R,K)
if ~isstruct(spin_system)||~isfield(spin_system,'model')||...
   ~isstruct(parameters)||~isscalar(parameters)||~isfield(parameters,'drive')
    error('finite_drive:input','Supply a retained subsystem and drive parameters.');
end
cce_pulse_events(parameters);
if isfield(parameters,'sequence')
    error('finite_drive:sequence','Use finite_drive native drive descriptors, not parameters.sequence.');
end
for k=1:numel(parameters.drive)
    ids=parameters.drive(k).spins;
    if ~isnumeric(ids)||~isreal(ids)||~isrow(ids)||isempty(ids)||...
       numel(unique(ids))~=numel(ids)||any(~ismember(ids,spin_system.model.observed_ids))
        error('finite_drive:spins','Finite drives address unique observed-spin identifiers.');
    end
end
if ~isnumeric(H)||~isnumeric(R)||~isnumeric(K)
    error('finite_drive:generators','Finite drives require conventional numeric generators.');
end
end


