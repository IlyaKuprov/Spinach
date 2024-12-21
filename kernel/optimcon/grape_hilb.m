function [traj_data,fidelity,grad] = grape_hilb(spin_system,drifts,controls,...
                                                waveform,rho_init,rho_targ,...
                                                fidelity_type) %#ok<*PFBNS>

% Hush up the output
shut_up.sys.output='hush';
shut_up.tols=spin_system.tols;
shut_up.sys.enable=spin_system.sys.enable;
shut_up.sys.disable=spin_system.sys.disable;

% Extract the timing grid
dt=spin_system.control.pulse_dt;

% Extract number of controls and time elements
nctrls = numel(controls);  nsteps = numel(dt);

% Make pointer arrays for trajectories
fwd_traj = cell(1,nsteps+1); fwd_traj{1}=rho_init;
bwd_traj = cell(1,nsteps+1); bwd_traj{1}=rho_targ;

% Preallocate gradient
grad = zeros([nctrls nsteps]);

% Reshape waveform
waveform=reshape(waveform,[nctrls nsteps]);

% Propagate
for n=1:nsteps

    % Decide current drifts
    if isscalar(drifts)

        % Time-independent drifts
        H_forw=drifts{1}; H_back=drifts{1};

    else

        % Time-dependent drifts
        H_forw=drifts{n}; H_back=drifts{nsteps+1-n};

    end

    % Add current controls to current drifts
    for k=1:nctrls

        % Forward evolution generator
        H_forw=H_forw+waveform(k,n)*controls{k};

        % Backward evolution generator
        H_back=H_back+waveform(k,nsteps+1-n)*controls{k};

    end

    % Propagators
    prop_forw=propagator(shut_up,H_forw,dt(n));
    prop_back=propagator(shut_up,H_back,-dt(nsteps+1-n));

    % Take the time step forward
    fwd_traj{n+1}=prop_forw*fwd_traj{n}*prop_forw';

    % Take the time step backward
    bwd_traj{n+1}=prop_back*bwd_traj{n}*prop_back';

end

% Compute fidelity
overlap = hdot(fwd_traj{end},rho_targ);

% Flip backwards trajectory
bwd_traj=fliplr(bwd_traj);

% Gradient loop
parfor n=1:nsteps

    % Decide current drifts
    if isscalar(drifts)

        % Time-independent drifts
        H_forw=drifts{1}; 

    else

        % Time-dependent drifts
        H_forw=drifts{n};

    end

    % Add current controls to current drifts
    for k=1:nctrls

        % Forward evolution generator
        H_forw=H_forw+waveform(k,n)*controls{k};

    end

    for k = 1:nctrls

        % Compute directional derivative w.r.t. control Matrix
        auxmat = dirdiff(spin_system,H_forw,controls{k},dt(n),2);

        % Compute Gradient
        grad(k,n) = real(hdot(bwd_traj{n+1},auxmat{2}*fwd_traj{n}*auxmat{1}'+...
                                            auxmat{1}*fwd_traj{n}*auxmat{2}'));

    end

end

% Fidelity and its derivatives
switch fidelity_type
    
    case {'real'}
        
        % Real part of the overlap
        fidelity=real(overlap);
        
        % Update gradient
        if exist('grad','var'), grad=real(grad); end
        
    case {'imag'}
        
        % Imaginary part of the overlap
        fidelity=imag(overlap);
        
        % Update gradient
        if exist('grad','var'), grad=imag(grad); end
        
    case {'square'}
        
        % Absolute square of the overalp
        fidelity=overlap*conj(overlap);
        
        % Product rule
        grad=grad*conj(overlap)+overlap*conj(grad);
        
        % Cleaning up
        grad=real(grad);
        
    otherwise
        
        % Complain and bomb out
        error('unknown fidelity type');
        
end

% Return trajectory data
traj_data.forward=fwd_traj;

end