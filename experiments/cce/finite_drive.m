% Finite rectangular or piecewise drive in a retained cluster.
% Syntax:
%
%     [answer,aux]=finite_drive(spin_system,parameters,H,R,K)
%
% parameters.times gives increasing elapsed seconds, zero or after the last
% finite slice ends. parameters.drive is a scalar structure with fields:
%   spins     - unique physical observed-spin identifiers
%   axis      - 'x', 'y', 'z', or a real unit three-vector
%   amplitude - real finite amplitude row in rad/s
%   detuning  - real finite longitudinal offset in rad/s
%   duration  - matching positive duration row in seconds
% Optional time specifies the start in seconds; otherwise the start is zero.
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
parameters.sequence=cce_pulse_events(parameters);
finite=strcmp({parameters.sequence.type},'finite');
stops=[parameters.sequence(finite).at]+[parameters.sequence(finite).duration];
if any(parameters.times>0&parameters.times<max(stops))
    error('finite_drive:times','Nonzero observations must follow completion of the finite drive.');
end

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
if isfield(parameters.drive,'spins')
    ids=parameters.drive.spins;
    if ~isnumeric(ids)||~isreal(ids)||~isrow(ids)||isempty(ids)||...
       numel(unique(ids))~=numel(ids)||any(~ismember(ids,spin_system.model.observed_ids))
        error('finite_drive:spins','Finite drives address unique observed-spin identifiers.');
    end
end
if ~isnumeric(H)||~isnumeric(R)||~isnumeric(K)
    error('finite_drive:generators','Finite drives require conventional numeric generators.');
end
end


