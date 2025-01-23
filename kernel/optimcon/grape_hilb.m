% Gradient Ascent Pulse Engineering (GRAPE) objective function and gradient
% Propagates the system through a user-supplied shaped pulse from a given 
% initial state and projects the result onto the given final  state. The 
% fidelity is returned, along with its gradient with respect to amplitudes 
% of all control operators at every time step of the shaped pulse. 
% Uses Hilbert-space formalism. Syntax:
%
%  [traj_data,fidelity,grad]=grape_hilb(spin_system,drifts,controls,...
%                                       waveform,rho_init,rho_targ,...
%                                       fidelity_type)
% Parameters:
%
%   spin_system         - Spinach data object that has been through 
%                         the optimcon.m problem setup function.
% 
%   drifts              - the drift Hamiltonians: a cell array con-
%                         taining one matrix (for time-independent 
%                         drift) or multiple matrices (for time-de-
%                         pendent drift).
%
%   controls            - control operators in Hilbert space (cell 
%                         array of matrices).
%
%   waveform            - control coefficients for each control ope-
%                         rator (in rows of a matrix), rad/s
%
%   rho_init            - initial state of the system as a matrix in
%                         Hilbert space.
%
%   rho_targ            - target state of the system as a matrix in
%                         Hilbert space.
%
%   fidelity_type       - 'real'   (real part of the overlap)
%                         'imag'   (imaginary part of the overlap)
%                         'square' (absolute square of the overlap)
%
% Outputs:
%
%   fidelity            - fidelity of the control sequence
%
%   grad                - gradient of the fidelity with respect to
%                         the control sequence
%
%   traj_data.forward   - forward trajectory from the initial condi-
%                         tion(a stack of state vectors)
%
% Note: this is a low level function that is not designed to be called 
%       directly. Use grape_xy.m and grape_phase.m instead.
%
% i.kuprov@soton.ac.uk
% m.keitel@soton.ac.uk
%
% TODO (Keitel): add logic to avoid computing backward trajectory 
%                when the gradient is not requested
%
% <https://spindynamics.org/wiki/index.php?title=grape.m>

function [traj_data,fidelity,grad] = grape_hilb(spin_system,drifts,controls,...
                                                waveform,rho_init,rho_targ,...
                                                fidelity_type) %#ok<*PFBNS>

% Check consistency
grumble(spin_system,drifts,controls,waveform,rho_init,rho_targ);

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
        H_forw=drifts{1}; 
        H_back=drifts{1};

    else

        % Time-dependent drifts
        H_forw=drifts{n}; 
        H_back=drifts{nsteps+1-n};

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

% Compute the Gradient
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

% Consistency enforcement
function grumble(spin_system,drifts,controls,waveform,rho_init,rho_targ)
if ~ismember(spin_system.bas.formalism,{'zeeman-hilb'})
    error('this function requires Lioville space formalism.');
end
if (~isnumeric(rho_init))||(~ismatrix(rho_init))||(size(rho_init,1)~=size(rho_init,2))
    error('rho_init must be a square matrix.');
end
if (~isnumeric(rho_targ))||(~ismatrix(rho_init))||(size(rho_targ,1)~=size(rho_targ,2))
    error('rho_targ must be a square vector.');
end
if ~iscell(drifts)
    error('drifts must be a cell array of matrices.');
end
for n=1:numel(drifts)     
    if (~isnumeric(drifts{n}))||(size(drifts{n},1)~=size(drifts{n},2))
        error('all elements of drifts cell array must be square matrices.');
    end
    if (size(drifts{n},1)~=size(rho_init,1))||...
       (size(drifts{n},1)~=size(rho_targ,1))
        error('dimensions of drift, rho_init and rho_targ must be consistent.');
    end
end
if ~iscell(controls)
    error('controls must be a cell array of square matrices.');
end
for n=1:numel(controls)
    if (~isnumeric(controls{n}))||...
       (size(controls{n},1)~=size(controls{n},2))||...
       (size(controls{n},1)~=size(drifts{1},1))
        error('control operators must have the same size as drift operators.');
    end
end
if (~isnumeric(waveform))||(~isreal(waveform))
    error('waveform must be a real numeric array.');
end
if size(waveform,1)~=numel(controls)
    error('number of waveform rows must be equal to the number of controls.');
end
if size(waveform,2)~=spin_system.control.pulse_ntpts
    error(['waveform must have ' int2str(spin_system.control.pulse_ntpts) ' columns.']);
end
end