% Computes a rigorous Hilbert-space reachability bound from Hamiltonian
% spectral-width limits. This is a deliberately conservative certificate
% intended for optimal-control benchmarks. Syntax:
%
%          certificate=reachable_bound(drift,controls,rho_init,rho_targ,...
%                                      amplitude_bounds,duration)
%
% Parameters:
%
%             drift  - Hilbert-space drift Hamiltonian, or a cell array of
%                      time-dependent drift Hamiltonians
%
%          controls  - cell array of Hilbert-space control Hamiltonians
%
%          rho_init  - initial state vector, state matrix, or observable
%
%          rho_targ  - target state vector, state matrix, or observable
%
%  amplitude_bounds  - control amplitude bounds, as a scalar absolute
%                      bound, a two-element [lower upper] vector, or a
%                      control-by-two array of [lower upper] rows
%
%          duration  - control horizon, seconds
%
% Outputs:
%
%     certificate.fidelity_upper_bound  - rigorous upper bound on the
%                                         target metric
%
%     certificate.time_lower_bound      - rigorous lower bound on the
%                                         time needed for the target
%
%     certificate.solver_status         - certificate construction status
%
%     certificate.relaxation_order      - certificate family identifier
%
%     certificate.angle_initial_target  - Fubini-Study or Hilbert-Schmidt
%                                         angle between source and target
%
%     certificate.angle_reachable       - largest certified angular
%                                         distance reachable in duration
%
%     certificate.overlap_upper_bound   - upper bound on normalised
%                                         target overlap
%
% For state vectors, the angle is the Fubini-Study angle and the fidelity
% bound is the squared state overlap. For matrices, the angle is the
% normalised Hilbert-Schmidt angle and the fidelity bound is the
% normalised real Hilbert-Schmidt overlap. The drift and bounded controls
% are combined by the triangle inequality, making the certificate
% rigorous but intentionally conservative. This function does not
% implement a moment/SOS hierarchy; the returned structure includes empty
% slots for future specialised SDP relaxations when tractable reductions
% are explicitly available.
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

% Compute a conservative Hamiltonian spectral-width ceiling
ham_width_bound=local_width_bound(drift);
for n=1:numel(controls)
    ham_width_bound=ham_width_bound+max(abs(amp_bounds(n,:)))*...
                    local_spectral_width(controls{n});
end

% Classify the source and target representation
target_type=local_target_type(rho_init,rho_targ);

% Compute the relevant Hilbert-space angle
switch target_type

    case 'state_vector'

        % Fubini-Study speed is bounded by half the spectral width
        source=rho_init(:)/norm(rho_init(:),2);
        target=rho_targ(:)/norm(rho_targ(:),2);
        init_overlap=abs(source'*target);
        angle_initial_target=acos(local_clip(init_overlap));
        angle_speed_bound=ham_width_bound/2;

    case 'operator'

        % Hilbert-Schmidt angular speed is bounded by the spectral width
        source_norm=sqrt(real(hdot(rho_init,rho_init)));
        target_norm=sqrt(real(hdot(rho_targ,rho_targ)));
        init_overlap=real(hdot(rho_init,rho_targ))/(source_norm*target_norm);
        angle_initial_target=acos(local_clip(init_overlap));
        angle_speed_bound=ham_width_bound;

    otherwise

        % Complain and bomb out
        error('unrecognised target representation.');

end

% Compute reachable angle and residual target gap
angle_reachable=angle_speed_bound*duration;
angle_gap=max(angle_initial_target-angle_reachable,0);
overlap_upper_bound=cos(angle_gap);

% Convert overlap bounds into the requested target metric
switch target_type

    case 'state_vector'

        % State-vector fidelity is the squared Fubini-Study overlap
        fidelity_upper_bound=overlap_upper_bound^2;

    case 'operator'

        % Matrix objective is the normalised real Hilbert-Schmidt overlap
        fidelity_upper_bound=overlap_upper_bound;

end

% Minimum time follows from the same angular-speed bound
if angle_speed_bound==0
    if angle_initial_target==0
        time_lower_bound=0;
    else
        time_lower_bound=Inf;
    end
else
    time_lower_bound=angle_initial_target/angle_speed_bound;
end

% Prepare slots for specialised reductions
sdp_relaxations=struct();
sdp_relaxations.single_excitation_chain=[];
sdp_relaxations.symmetry_reduced_few_spin=[];

% Return a certificate structure
certificate=struct();
certificate.fidelity_upper_bound=fidelity_upper_bound;
certificate.time_lower_bound=time_lower_bound;
certificate.solver_status='closed_form_speed_limit';
certificate.relaxation_order=0;
certificate.angle_initial_target=angle_initial_target;
certificate.angle_reachable=angle_reachable;
certificate.overlap_upper_bound=overlap_upper_bound;
certificate.hamiltonian_spectral_width_bound=ham_width_bound;
certificate.angle_speed_bound=angle_speed_bound;
certificate.duration=duration;
certificate.amplitude_bounds=amp_bounds;
certificate.target_type=target_type;
certificate.scope='arbitrary-dimension Hilbert-space speed-limit certificate';
certificate.sdp_relaxation_status='not_attempted';
certificate.specialised_sdp_relaxations=sdp_relaxations;

end

% Consistency enforcement
function grumble(drift,controls,rho_init,rho_targ,amplitude_bounds,duration)
matrix_dim=local_drift_dim(drift);
local_check_drift(drift,matrix_dim);
if ~iscell(controls)
    error('controls must be a cell array of numeric matrices.');
end
for n=1:numel(controls)
    if (~isnumeric(controls{n}))||(~isequal(size(controls{n}),[matrix_dim matrix_dim]))
        error('all controls must be numeric matrices matching drift dimensions.');
    end
    if any(~isfinite(controls{n}(:)))
        error('all controls must contain finite numbers.');
    end
    if cheap_norm(controls{n}-controls{n}')>1e-10*max(cheap_norm(controls{n}),1)
        error('all controls must be Hermitian.');
    end
end
if (~isnumeric(rho_init))||(~isnumeric(rho_targ))
    error('rho_init and rho_targ must be numeric.');
end
if any(~isfinite(rho_init(:)))||any(~isfinite(rho_targ(:)))
    error('rho_init and rho_targ must contain finite numbers.');
end
if ~local_valid_target(rho_init,rho_targ,matrix_dim)
    error('rho_init and rho_targ must be matching state vectors or square matrices.');
end
if local_target_norm(rho_init)==0
    error('rho_init must have a non-zero norm.');
end
if local_target_norm(rho_targ)==0
    error('rho_targ must have a non-zero norm.');
end
if (~isnumeric(amplitude_bounds))||(~isreal(amplitude_bounds))||...
   isempty(amplitude_bounds)||any(~isfinite(amplitude_bounds(:)))
    error('amplitude_bounds must be a non-empty finite real numeric array.');
end
if (~isnumeric(duration))||(~isreal(duration))||(~isscalar(duration))||...
   (~isfinite(duration))||(duration<0)
    error('duration must be a non-negative finite real scalar.');
end
end

% Return Hilbert-space dimension from drift storage
function matrix_dim=local_drift_dim(drift)
if iscell(drift)
    if isempty(drift)
        error('drift must not be an empty cell array.');
    end
    matrix_dim=local_drift_dim(drift{1});
else
    if (~isnumeric(drift))||(~ismatrix(drift))||(size(drift,1)~=size(drift,2))
        error('drift must be a numeric square matrix or a cell array thereof.');
    end
    matrix_dim=size(drift,1);
end
end

% Validate drift matrices in nested storage
function local_check_drift(drift,matrix_dim)
if iscell(drift)
    if isempty(drift)
        error('drift must not be an empty cell array.');
    end
    for n=1:numel(drift)
        local_check_drift(drift{n},matrix_dim);
    end
else
    if (~isnumeric(drift))||(~isequal(size(drift),[matrix_dim matrix_dim]))
        error('all drift matrices must be numeric matrices of matching dimensions.');
    end
    if any(~isfinite(drift(:)))
        error('all drift matrices must contain finite numbers.');
    end
    if cheap_norm(drift-drift')>1e-10*max(cheap_norm(drift),1)
        error('all drift matrices must be Hermitian.');
    end
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

% Conservative spectral-width bound for nested drift storage
function width_bound=local_width_bound(drift)
if iscell(drift)
    width_bound=0;
    for n=1:numel(drift)
        width_bound=max(width_bound,local_width_bound(drift{n}));
    end
else
    width_bound=local_spectral_width(drift);
end
end

% Spectral width of a Hermitian matrix
function width=local_spectral_width(H)
lambda=eig(full((H+H')/2));
width=max(real(lambda))-min(real(lambda));
end

% Classify source and target representation
function target_type=local_target_type(rho_init,rho_targ)
if isvector(rho_init)&&isvector(rho_targ)
    target_type='state_vector';
else
    target_type='operator';
end
end

% Validate source and target dimensions
function verdict=local_valid_target(rho_init,rho_targ,matrix_dim)
if isvector(rho_init)&&isvector(rho_targ)
    verdict=(numel(rho_init)==matrix_dim)&&(numel(rho_targ)==matrix_dim);
else
    verdict=isequal(size(rho_init),[matrix_dim matrix_dim])&&...
            isequal(size(rho_targ),[matrix_dim matrix_dim]);
end
end

% Compute vector or matrix norm
function rho_norm=local_target_norm(rho)
if isvector(rho)
    rho_norm=norm(rho(:),2);
else
    rho_norm=sqrt(real(hdot(rho,rho)));
end
end

% Clip overlap into the arccos domain
function value=local_clip(value)
value=max(min(real(value),1),-1);
end
