% Certifies a small Spinach optimal-control problem against a rigorous
% reachability bound. Syntax:
%
%             certificate=oc_certify(spin_system,problem)
%
% Parameters:
%
%     spin_system  - Spinach data object that has been through
%                    optimcon.m, or a minimal object with control
%                    fields populated by the caller
%
%        problem   - optional structure overriding fields extracted
%                    from spin_system.control:
%
%                      drift             2x2 drift Hamiltonian
%                      controls          cell array of 2x2 controls
%                      rho_init          initial state or observable
%                      rho_targ          target state or observable
%                      amplitude_bounds  scalar, [lower upper], or
%                                        control-by-two bounds
%                      duration          total pulse duration
%                      waveform          optional GRAPE waveform
%                      fidelity_type     GRAPE fidelity type
%
% Outputs:
%
%     certificate  - structure returned by reachable_bound(), augmented
%                    with best_fidelity and fidelity_gap fields when a
%                    waveform is supplied
%
% This is an experimental certification front end. The first supported
% certificate is intentionally narrow: two-level Hilbert-space systems
% with bounded Hamiltonian amplitudes. It is useful as a GRAPE benchmark
% for one-spin examples and as a scaffold for later SDP-backed moment
% relaxations.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=oc_certify.m>

function certificate=oc_certify(spin_system,problem)

% Default problem structure
if nargin<2, problem=struct(); end

% Check consistency
grumble(spin_system,problem);

% Pull problem data from spin_system.control unless overridden
control=spin_system.control;
if isfield(problem,'drift')
    drift=problem.drift;
else
    drift=local_first_drift(control.drifts);
end
if isfield(problem,'controls')
    controls=problem.controls;
else
    controls=control.operators;
end
if isfield(problem,'rho_init')
    rho_init=problem.rho_init;
else
    rho_init=control.rho_init{1};
end
if isfield(problem,'rho_targ')
    rho_targ=problem.rho_targ;
else
    rho_targ=control.rho_targ{1};
end
if isfield(problem,'duration')
    duration=problem.duration;
else
    duration=sum(control.pulse_dt);
end
if isfield(problem,'amplitude_bounds')
    amplitude_bounds=problem.amplitude_bounds;
else
    amplitude_bounds=local_control_bounds(control,numel(controls));
end

% Get the reachability certificate
certificate=reachable_bound(drift,controls,rho_init,rho_targ,...
                            amplitude_bounds,duration);

% Add optional GRAPE comparison
certificate.best_fidelity=[];
certificate.fidelity_gap=[];
if isfield(problem,'waveform')
    if isfield(problem,'fidelity_type')
        fidelity_type=problem.fidelity_type;
    elseif isfield(control,'fidelity')
        fidelity_type=control.fidelity;
    else
        fidelity_type='real';
    end
    [~,best_fidelity]=grape_hilb(spin_system,{drift},controls,...
                                 problem.waveform,rho_init,rho_targ,...
                                 fidelity_type);
    certificate.best_fidelity=best_fidelity;
    certificate.fidelity_gap=certificate.overlap_upper_bound-best_fidelity;
end

end

% Consistency enforcement
function grumble(spin_system,problem)
if (~isstruct(spin_system))||(~isfield(spin_system,'control'))
    error('spin_system must be a Spinach object with a control field.');
end
if ~isstruct(problem)
    error('problem must be a structure.');
end
end

% Extract first drift matrix from Spinach's nested control.drifts layout
function drift=local_first_drift(drifts)
if iscell(drifts)
    drift=drifts{1};
    while iscell(drift)
        drift=drift{1};
    end
else
    drift=drifts;
end
end

% Infer amplitude bounds from optimcon control metadata
function amplitude_bounds=local_control_bounds(control,n_controls)
if isfield(control,'l_bound')&&isfield(control,'u_bound')
    amplitude_bounds=repmat([control.l_bound control.u_bound],n_controls,1);
elseif isfield(control,'pwr_levels')
    amplitude_bounds=max(abs(control.pwr_levels(:)));
else
    error('amplitude bounds are required: supply problem.amplitude_bounds or control bounds.');
end
end

