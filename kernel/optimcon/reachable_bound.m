% Computes a rigorous two-level reachability bound from Hamiltonian
% norm limits. This is a deliberately small certificate intended for
% pedagogical optimal-control benchmarks. Syntax:
%
%          certificate=reachable_bound(drift,controls,rho_init,rho_targ,...
%                                      amplitude_bounds,duration)
%
% Parameters:
%
%             drift  - 2x2 Hilbert-space drift Hamiltonian
%
%          controls  - cell array of 2x2 Hilbert-space control
%                      Hamiltonians
%
%          rho_init  - initial two-level state or observable
%
%          rho_targ  - target two-level state or observable
%
%  amplitude_bounds  - control amplitude bounds, as a scalar absolute
%                      bound, a two-element [lower upper] vector, or
%                      a control-by-two array of [lower upper] rows
%
%          duration  - control horizon, seconds
%
% Outputs:
%
%     certificate.fidelity_upper_bound  - rigorous upper bound on the
%                                         reachable target fidelity
%
%     certificate.time_lower_bound      - rigorous lower bound on the
%                                         time needed for the target
%
%     certificate.solver_status         - certificate construction status
%
%     certificate.relaxation_order      - certificate family identifier
%
%     certificate.angle_initial_target  - Bloch-sphere angle between
%                                         source and target
%
%     certificate.angle_reachable       - largest certified angular
%                                         distance reachable in duration
%
%     certificate.overlap_upper_bound   - upper bound on normalised
%                                         Hilbert-Schmidt overlap
%
% This function does not implement a general moment/SOS hierarchy. It
% supplies a small, closed-form certificate for one-qubit and one-spin
% Hilbert-space benchmarks: the angular velocity on the Bloch sphere is
% bounded by the spectral width of the Hamiltonian. Drift and bounded
% controls are combined by the triangle inequality, making the resulting
% bound rigorous but intentionally conservative.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=reachable_bound.m>

function certificate=reachable_bound(drift,controls,rho_init,rho_targ,...
                                     amplitude_bounds,duration)

% Check consistency
grumble(drift,controls,rho_init,rho_targ,amplitude_bounds,duration);

% Put amplitude limits into one row per control
amp_bounds=local_amplitude_bounds(amplitude_bounds,numel(controls));

% Compute a rigorous angular-speed ceiling
omega_bound=local_spectral_width(drift);
for n=1:numel(controls)
    omega_bound=omega_bound+max(abs(amp_bounds(n,:)))*...
                local_spectral_width(controls{n});
end

% Convert source and target to Bloch vectors
source=local_bloch_vector(rho_init);
target=local_bloch_vector(rho_targ);

% Normalise the Bloch vectors
source_norm=norm(source,2);
target_norm=norm(target,2);
if (source_norm==0)||(target_norm==0)
    error('rho_init and rho_targ must have non-zero Bloch vectors.');
end
source=source/source_norm;
target=target/target_norm;

% Initial target angle and the maximum angle that can be swept
angle_initial_target=acos(max(min(real(source'*target),1),-1));
angle_reachable=omega_bound*duration;
angle_gap=max(angle_initial_target-angle_reachable,0);

% Normalised overlap and pure-state fidelity upper bounds
overlap_upper_bound=cos(angle_gap);
fidelity_upper_bound=(1+overlap_upper_bound)/2;

% Minimum time follows from the same angular-speed bound
if omega_bound==0
    if angle_initial_target==0
        time_lower_bound=0;
    else
        time_lower_bound=Inf;
    end
else
    time_lower_bound=angle_initial_target/omega_bound;
end

% Return a certificate structure
certificate=struct();
certificate.fidelity_upper_bound=fidelity_upper_bound;
certificate.time_lower_bound=time_lower_bound;
certificate.solver_status='closed_form_speed_limit';
certificate.relaxation_order=0;
certificate.angle_initial_target=angle_initial_target;
certificate.angle_reachable=angle_reachable;
certificate.overlap_upper_bound=overlap_upper_bound;
certificate.hamiltonian_spectral_width_bound=omega_bound;
certificate.duration=duration;
certificate.amplitude_bounds=amp_bounds;
certificate.scope='two-level Hilbert-space speed-limit certificate';

end

% Consistency enforcement
function grumble(drift,controls,rho_init,rho_targ,amplitude_bounds,duration)
if (~isnumeric(drift))||(~isequal(size(drift),[2 2]))
    error('drift must be a 2x2 numeric matrix.');
end
if cheap_norm(drift-drift')>1e-10*max(cheap_norm(drift),1)
    error('drift must be Hermitian.');
end
if ~iscell(controls)
    error('controls must be a cell array of 2x2 matrices.');
end
for n=1:numel(controls)
    if (~isnumeric(controls{n}))||(~isequal(size(controls{n}),[2 2]))
        error('all controls must be 2x2 numeric matrices.');
    end
    if cheap_norm(controls{n}-controls{n}')>1e-10*max(cheap_norm(controls{n}),1)
        error('all controls must be Hermitian.');
    end
end
if (~isnumeric(rho_init))||(~isequal(size(rho_init),[2 2]))||...
   (~isnumeric(rho_targ))||(~isequal(size(rho_targ),[2 2]))
    error('rho_init and rho_targ must be 2x2 numeric matrices.');
end
if (~isnumeric(amplitude_bounds))||(~isreal(amplitude_bounds))||isempty(amplitude_bounds)
    error('amplitude_bounds must be a non-empty real numeric array.');
end
if (~isnumeric(duration))||(~isreal(duration))||(~isscalar(duration))||(duration<0)
    error('duration must be a non-negative real scalar.');
end
end

% Normalise amplitude bound formats
function amp_bounds=local_amplitude_bounds(amplitude_bounds,n_controls)
if isscalar(amplitude_bounds)
    amp_bounds=repmat([-abs(amplitude_bounds) abs(amplitude_bounds)],n_controls,1);
elseif isvector(amplitude_bounds)&&(numel(amplitude_bounds)==2)
    amp_bounds=repmat(amplitude_bounds(:).',n_controls,1);
elseif isequal(size(amplitude_bounds),[n_controls 2])
    amp_bounds=amplitude_bounds;
else
    error('amplitude_bounds must be scalar, [lower upper], or control-by-two.');
end
if any(amp_bounds(:,1)>amp_bounds(:,2))
    error('each lower amplitude bound must not exceed the upper bound.');
end
end

% Spectral width of a Hermitian matrix
function width=local_spectral_width(H)
lambda=eig(full((H+H')/2));
width=max(real(lambda))-min(real(lambda));
end

% Bloch vector of a two-level operator
function bloch=local_bloch_vector(rho)
S=pauli(2);
sigma={2*S.x,2*S.y,2*S.z};
bloch=zeros(3,1);
rho=(rho+rho')/2;
for n=1:3
    bloch(n)=real(trace(rho*sigma{n}));
end
end

