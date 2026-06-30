% Certifies a Spinach Hilbert-space optimal-control problem against a
% rigorous reachability bound. Syntax:
%
%             certificate=oc_certify(spin_system,problem)
%
% Parameters:
%
%     spin_system  - Spinach data object that has been through
%                    optimcon.m, or a minimal object with control fields
%                    populated by the caller
%
%        problem   - structure overriding fields extracted from
%                    spin_system.control:
%
%                      drift             Hilbert-space drift Hamiltonian,
%                                        or cell array of time slices
%                      controls          cell array of control Hamiltonians
%                      rho_init          initial state or observable
%                      rho_targ          target state or observable
%                      amplitude_bounds  scalar, [lower upper], or
%                                        control-by-two bounds
%                      duration          total pulse duration
%                      waveform          candidate optimcon waveform
%
% Outputs:
%
%     certificate  - structure returned by reachable_bound(), augmented
%                    with best_fidelity and fidelity_gap fields when a
%                    waveform is supplied
%
% This is an experimental certification front end. The first supported
% certificate is intentionally conservative: arbitrary-dimensional
% Hilbert-space systems with bounded Hamiltonian amplitudes. It gives a
% closed-form speed-limit certificate and reserves structure slots for
% future specialised SDP relaxations when a tractable reduction is
% explicitly available.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=oc_certify.m>

function certificate=oc_certify(spin_system,problem)

% Check consistency
grumble(spin_system,problem);

% Pull problem data from spin_system.control unless overridden
control=spin_system.control;
if isfield(problem,'drift')
    drift=problem.drift;
else
    drift=local_control_drifts(control);
end
if isfield(problem,'controls')
    controls=problem.controls;
else
    controls=control.operators;
end
if xor(isfield(problem,'rho_init'),isfield(problem,'rho_targ'))
    error('problem.rho_init and problem.rho_targ must be supplied together.');
end
if isfield(problem,'rho_init')
    rho_init=problem.rho_init;
    rho_targ=problem.rho_targ;
else
    if (numel(control.rho_init)~=1)||(numel(control.rho_targ)~=1)
        error('oc_certify currently requires one state-target pair.');
    end
    rho_init=control.rho_init{1};
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

% Add optional waveform comparison on the same normalised scale
certificate.best_fidelity=[];
certificate.fidelity_gap=[];
certificate.best_overlap=[];
certificate.overlap_gap=[];
if isfield(problem,'waveform')
    waveform=local_physical_waveform(problem.waveform,control);
    [best_fidelity,best_overlap]=local_waveform_metric(spin_system,drift,...
                                                       controls,waveform,...
                                                       rho_init,rho_targ,...
                                                       certificate.target_type);
    certificate.best_fidelity=best_fidelity;
    certificate.fidelity_gap=certificate.fidelity_upper_bound-best_fidelity;
    certificate.best_overlap=best_overlap;
    certificate.overlap_gap=certificate.overlap_upper_bound-best_overlap;
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

% Extract drift data from control metadata or value-store labels
function drift=local_control_drifts(control)
if isfield(control,'drifts')
    drift=control.drifts;
elseif isfield(control,'ndrifts')
    pool=gcp('nocreate');
    if isempty(pool)
        error('control.drifts are unavailable; supply problem.drift explicitly.');
    end
    store=pool.ValueStore;
    drift=cell(1,control.ndrifts);
    for n=1:control.ndrifts
        drift{n}=store(['oc_drift_' num2str(n)]);
    end
else
    error('control.drifts or problem.drift must be supplied.');
end
end

% Infer physical amplitude bounds from optimcon control metadata
function amplitude_bounds=local_control_bounds(control,n_controls)
if isfield(control,'pwr_levels')
    power_bound=max(abs(control.pwr_levels(:)));
else
    power_bound=1;
end
if isfield(control,'l_bound')&&isfield(control,'u_bound')
    amplitude_bounds=zeros(n_controls,2);
    for n=1:n_controls
        lower_slice=local_bound_slice(control.l_bound,n,n_controls);
        upper_slice=local_bound_slice(control.u_bound,n,n_controls);
        amp_bound=power_bound*max(abs([lower_slice(:); upper_slice(:)]));
        amplitude_bounds(n,:)=[-amp_bound amp_bound];
    end
elseif isfield(control,'pwr_levels')
    amplitude_bounds=repmat([-power_bound power_bound],n_controls,1);
else
    error('amplitude bounds are required: supply problem.amplitude_bounds or control bounds.');
end
end

% Select one control row from scalar or waveform-shaped bounds
function bound_slice=local_bound_slice(bounds,control_idx,n_controls)
if isscalar(bounds)
    bound_slice=bounds;
elseif size(bounds,1)==n_controls
    bound_slice=bounds(control_idx,:);
else
    error('control bounds must be scalar or have one row per control.');
end
end

% Convert an optimcon waveform into physical units
function waveform=local_physical_waveform(waveform,control)
if (~isnumeric(waveform))||(~isreal(waveform))||any(~isfinite(waveform(:)))
    error('problem.waveform must be a finite real numeric array.');
end
if isfield(control,'pwr_levels')
    if numel(control.pwr_levels)~=1
        error('waveform comparison currently requires one control power level.');
    end
    waveform=control.pwr_levels(1)*waveform;
end
end

% Evaluate a waveform without GRAPE optimiser zero-overlap guards
function [best_fidelity,best_overlap]=local_waveform_metric(spin_system,drift,...
                                                            controls,waveform,...
                                                            rho_init,rho_targ,...
                                                            target_type)

% Check waveform comparison metadata
control=spin_system.control;
if ~isfield(control,'pulse_dt')
    error('spin_system.control.pulse_dt is required for waveform comparison.');
end
if ~isfield(control,'integrator')
    error('spin_system.control.integrator is required for waveform comparison.');
end
if (size(waveform,1)~=numel(controls))
    error('problem.waveform must have one row per control operator.');
end

% Reject unsupported trajectory modifiers
if isfield(control,'dead_time')&&(control.dead_time~=0)
    error('waveform comparison does not currently support dead time.');
end
if isfield(control,'prefix')&&(~isempty(control.prefix))
    error('waveform comparison does not currently support prefix pulses.');
end
if isfield(control,'suffix')&&(~isempty(control.suffix))
    error('waveform comparison does not currently support suffix pulses.');
end

% Extract the drift time series
drift_series=local_drift_series(drift);
dt=control.pulse_dt;
rho=rho_init;

% Propagate with the requested time integrator
switch control.integrator

    case 'rectangle'

        % Check rectangle waveform dimensions
        if size(waveform,2)~=numel(dt)
            error('problem.waveform must have one column per time interval.');
        end

        % Run piecewise-constant propagation
        for n=1:numel(dt)
            H=drift_series{mod(n-1,numel(drift_series))+1};
            for k=1:numel(controls)
                H=H+waveform(k,n)*controls{k};
            end
            P=propagator(spin_system,H,dt(n));
            rho=local_apply_prop(P,rho,target_type);
            rho=local_apply_keyhole(control,rho,n);
        end

    case 'trapezium'

        % Check trapezium waveform dimensions
        if size(waveform,2)~=numel(dt)+1
            error('problem.waveform must have one more column than the time grid.');
        end

        % Run piecewise-linear propagation
        for n=1:numel(dt)
            left_ham=drift_series{mod(n-1,numel(drift_series))+1};
            right_ham=drift_series{mod(n,numel(drift_series))+1};
            for k=1:numel(controls)
                left_ham=left_ham+waveform(k,n)*controls{k};
                right_ham=right_ham+waveform(k,n+1)*controls{k};
            end
            H=(left_ham+right_ham)/2+...
              1i*dt(n)*(1/12)*(left_ham*right_ham-right_ham*left_ham);
            P=propagator(spin_system,H,dt(n));
            rho=local_apply_prop(P,rho,target_type);
            rho=local_apply_keyhole(control,rho,n);
        end

    otherwise

        % Complain and bomb out
        error('unknown integrator: must be ''rectangle'' or ''trapezium''.');

end

% Return the normalised target metric
[best_fidelity,best_overlap]=local_target_metric(rho,rho_targ,target_type);

end

% Extract one drift time series for waveform comparison
function drift_series=local_drift_series(drift)
if isnumeric(drift)
    drift_series={drift};
elseif all(cellfun(@isnumeric,drift(:)))
    drift_series=drift(:).';
elseif isscalar(drift)
    drift_series=local_drift_series(drift{1});
else
    error('waveform comparison currently requires one drift ensemble member.');
end
end

% Apply Hilbert-space propagation to vectors or operators
function rho=local_apply_prop(P,rho,target_type)
switch target_type

    case 'state_vector'

        % Propagate a state vector
        rho=P*rho;

    case 'operator'

        % Propagate a density matrix or observable
        rho=P*rho*P';

end
end

% Apply a keyhole function when present
function rho=local_apply_keyhole(control,rho,time_idx)
if isfield(control,'keyholes')&&(numel(control.keyholes)>=time_idx)&&...
   (~isempty(control.keyholes{time_idx}))
    rho=control.keyholes{time_idx}(rho);
end
end

% Compute the normalised target metric
function [best_fidelity,best_overlap]=local_target_metric(rho,rho_targ,target_type)
switch target_type

    case 'state_vector'

        % State-vector fidelity is the squared Fubini-Study overlap
        best_overlap=abs(rho_targ(:)'*rho(:))/(norm(rho(:),2)*norm(rho_targ(:),2));
        best_overlap=max(min(real(best_overlap),1),0);
        best_fidelity=best_overlap^2;

    case 'operator'

        % Matrix objective is the normalised real Hilbert-Schmidt overlap
        rho_norm=sqrt(real(hdot(rho,rho)));
        targ_norm=sqrt(real(hdot(rho_targ,rho_targ)));
        best_overlap=real(hdot(rho,rho_targ))/(rho_norm*targ_norm);
        best_overlap=max(min(real(best_overlap),1),-1);
        best_fidelity=best_overlap;

end
end
