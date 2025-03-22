% Gradient Ascent Pulse Engineering (GRAPE) objective function and gradient
% Propagates the system through a user-supplied shaped pulse from a given 
% initial state and projects the result onto the given final state. The fi-
% delity is returned, along with its gradient with respect to amplitudes of
% all control operators at every time step of the shaped pulse. This func-
% tion is for the Hilbert space version of GRAPE. Syntax:
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
%   rho_init            - initial state of the system as a density 
%                         matrix in Hilbert space.
%
%   rho_targ            - target state of the system as a density
%                         matrix in Hilbert space.
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
% ilya.kuprov@weizmann.ac.il
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=grape_hilb.m>

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
if nargout>2
    bwd_traj = cell(1,nsteps+1); bwd_traj{1}=rho_targ;
end

% Preallocate gradient
if nargout>2
    grad = zeros([nctrls nsteps]);
end

% Reshape waveform
waveform=reshape(waveform,[nctrls nsteps]);

% Make generators and propagators
H=cell(nsteps,1); P=cell(nsteps,1);
parfor n=1:nsteps

    % Decide current drifts
    if isscalar(drifts)

        % Time-independent drifts
        H{n}=drifts{1};
      
    else

        % Time-dependent drifts
        H{n}=drifts{n};
      
    end

    % Add current controls
    for k=1:nctrls

        % Forward evolution generator
        H{n}=H{n}+waveform(k,n)*controls{k};

    end

    % Make sure generator is Hermitian
    if cheap_norm(H{n}-H{n}')>1e-6

        % Bomb out if significantly non-Hermitian
        error('evolution generator must be Hermitian.');

    end

    % Tidy up generator and compute propagator
    P{n}=propagator(shut_up,(H{n}+H{n}')/2,dt(n));

end

% Forward trajectory
for n=1:nsteps
    fwd_traj{n+1}=P{n}*fwd_traj{n}*P{n}';
end

% Backward trajectory
if nargout>2
    for n=1:nsteps
        bwd_traj{n+1}=P{nsteps+1-n}'*bwd_traj{n}*P{nsteps+1-n};
    end
    bwd_traj=fliplr(bwd_traj);
end

% Compute fidelity
overlap=hdot(fwd_traj{end},rho_targ);

% Gradient loop
if nargout>2
    parfor n=1:nsteps
        for k=1:nctrls

            % Compute directional derivative w.r.t. control
            auxmat=dirdiff(shut_up,H{n},controls{k},dt(n),2);

            % Compute gradient element
            grad(k,n)=hdot(bwd_traj{n+1},auxmat{2}*fwd_traj{n}*auxmat{1}'+...
                                         auxmat{1}*fwd_traj{n}*auxmat{2}');

        end
    end
end

% Fidelity and its derivatives
switch fidelity_type
    
    case {'real'}
        
        % Real part of the overlap
        fidelity=real(overlap);
        
        % Update gradient
        if nargout>2, grad=real(grad); end
        
    case {'imag'}
        
        % Imaginary part of the overlap
        fidelity=imag(overlap);
        
        % Update gradient
        if nargout>2, grad=imag(grad); end
        
    case {'square'}
        
        % Absolute square of the overalp
        fidelity=overlap*conj(overlap);
        
        % Update gradient
        if nargout>2
            grad=real(grad*conj(overlap)+...
                      overlap*conj(grad));
        end  
      
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
    error('this function requires a density matrix based formalism.');
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

% Series of headlines in Le Moniteur Universel ("le journal de la pensée 
% officielle"), as Napoleon was approaching Paris between 9 and 21 March
% 1815, according to Alexandre Dumas:
%
%  "L'anthropophage est sorti de son repaire."
%
%  "L'ogre de Corse vient de débarquer au golfe Juan."
%
%  "Le tigre est arrivé à Gap."
%
%  "Le monstre a couché à Grenoble."
%
%  "Le tyran a traversé Lyon."
%
%  "L'usurpateur a été vu à soixante lieues de la capitale."
%
%  "Bonaparte s'avance à grands pas, mais il n'entrera jamais dans Paris."
%
%  "Napoléon sera demain sous nos remparts."
%
%  "L'empereur est arrivé à Fontainebleau."
%
%  "Sa Majesté Impériale et Royale a fait hier son entrée en son château 
%   des Tuileries au milieu de ses fidèles sujets."
%
% Dumas concludes: "This is the ultimate monument to journalism. It need not
% do anything else, for it won't do anything better."

