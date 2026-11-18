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
% This experiment has no controls. Use finite_drive for finite pulses,
% refocusing_train for ideal trains, or hahn_echo for a Hahn experiment.
% Caller-supplied sequence fields are rejected even when empty.
% H,R,K are conventional angular-frequency generators or tagged projected
% branches from cce_generators. L=H+1i*R+1i*K.
%
% ilya.kuprov@weizmann.ac.il

function [answer,aux]=free_evolution(spin_system,parameters,H,R,K)

% Validate this experiment before shared preparation, evolution, and detection
grumble(spin_system,parameters,H,R,K);
[answer,aux]=cce_acquire(spin_system,parameters,H,R,K);

end

% Enforce free evolution without caller-supplied pulse specifications
function grumble(spin_system,parameters,H,R,K)
if ~isstruct(spin_system)||~isstruct(parameters)||~isscalar(parameters)||...
   ~(isnumeric(H)||isstruct(H))||~isnumeric(R)||~isnumeric(K)
    error('free_evolution:input','Supply a native subsystem, parameters, and generators.');
end
if any(isfield(parameters,{'sequence','drive','refocusing','train'}))
    error('free_evolution:sequence','Use finite_drive, refocusing_train, or hahn_echo for controls.');
end
end


