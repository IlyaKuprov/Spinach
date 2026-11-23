% Free evolution of a retained cluster with experiment-owned detection.
% Syntax:
%
%     [answer,aux]=free_evolution(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    spin_system - native retained subsystem from cce_subsystem
%
%    parameters  - experiment structure with the fields described below
%
%    H,R,K       - angular-frequency generators in the retained basis;
%                  L=H+1i*R+1i*K for conventional matrix generators
%
% parameters.times gives elapsed times in seconds. rho0 and coil are explicit
% observed-spin matrices, or transition selects preparation and detection.
% preparation may contain state and detector descriptors. normalized defaults
% to true for signals. output is signal, observable, or density; density output
% is the column-vectorised observed density, one row per elapsed time.
% Optional aux retains coherent spectral residuals and absolute error bounds.
% This experiment has no controls. Use finite_drive for finite pulses,
% refocusing_train for ideal trains, or hahn_echo for a Hahn experiment.
% Caller-supplied sequence fields are rejected even when empty.
% H,R,K are conventional angular-frequency generators, or the projected
% branch pair and dissipator from cce_generators. L=H+1i*R+1i*K.
%
% Outputs:
%
%    answer - one row per elapsed time: detected signals or observables,
%             or the column-vectorised observed density
%
%    aux    - propagation diagnostics and absolute arithmetic error bounds
%
% ilya.kuprov@weizmann.ac.il

function [answer,aux]=free_evolution(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Prepare, propagate, and detect the retained subsystem
[answer,aux]=cce_acquire(spin_system,parameters,H,R,K);

end

% Enforce free evolution without caller-supplied pulse specifications
function grumble(spin_system,parameters,H,R,K)
if ~isstruct(parameters)||~isscalar(parameters)||~isfield(parameters,'times')||...
   ~isnumeric(parameters.times)||~isreal(parameters.times)||~isvector(parameters.times)||...
   any(~isfinite(parameters.times))||any(parameters.times<0)||any(diff(parameters.times)<0)
    error('free_evolution:times','Supply finite nonnegative increasing elapsed times.');
end
if ~isstruct(spin_system)||~isstruct(parameters)||~isscalar(parameters)||...
   ~(isnumeric(H)||iscell(H))||~isnumeric(R)||~isnumeric(K)
    error('free_evolution:input',...
          'Supply a native subsystem, parameters, and generators.');
end
if any(isfield(parameters,{'sequence','drive','refocusing','train'}))
    error('free_evolution:sequence',...
          'Use finite_drive, refocusing_train, or hahn_echo for controls.');
end
end


