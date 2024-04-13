% A special case of Gradient Ascent Pulse Engineering (GRAPE) objective
% function and gradient with respect to the vector of waveform slice du-
% rations. Syntax:
%
%    [fidelity,grad]=tgrape(spin_system,drift,controls,waveform,...
%                           dt_grid,time_unit,rho_init,rho_targ)
%
% Parameters:
%
%   drift               - the drift Liouvillian, a matrix
%
%   controls            - control operators, a cell array 
%                         of matrices
%
%   waveform            - control coefficients for each control 
%                         operator (columns) at each time slice
%                         (rows), rad/s
%
%   dt_grid             - time slice durations, a col vector
%                         in the units of time chosen so that
%                         the elements are of the order of 1
%
%   time_unit           - unit of time, seconds; this is needed
%                         because optimisers get stuck when the
%                         variables are badly scaled
%   
%   rho_init            - initial state of the system, a column
%                         vector
%
%   rho_targ            - target state of the system, a column
%                         vector
%
% Outputs:
%
%   fidelity            - fidelity of the control sequence
%
%   grad                - gradient of the fidelity with respect to
%                         the durations of the slices
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=tgrape.m>

function [fidelity,grad]=tgrape(spin_system,drift,controls,waveform,...
                                dt_grid,time_unit,rho_init,rho_targ)
% Check consistency
grumble(spin_system,drift,controls,waveform,dt_grid,rho_init,rho_targ);
    
% Count the time steps
nsteps=size(waveform,2);

% Convert time units
dt_grid=dt_grid*time_unit;

% Preallocate trajectories
fwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i);
bwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i);

% Hush up the output
spin_system.sys.output='hush';

% Initialise forward and backward trajectories
fwd_traj(:,1)=rho_init; bwd_traj(:,1)=rho_targ;

% Compute trajectories
for n=1:nsteps

    % Start with the drift generator
    L_forw=drift; L_back=drift';

    % Add current controls
    for k=1:numel(controls)
        L_forw=L_forw+waveform(k,n)*controls{k};
        L_back=L_back+waveform(k,nsteps+1-n)*controls{k};
    end

    % Take time steps forwards and backwards
    fwd_traj(:,n+1)=step(spin_system,L_forw,fwd_traj(:,n),+dt_grid(n));
    bwd_traj(:,n+1)=step(spin_system,L_back,bwd_traj(:,n),-dt_grid(nsteps+1-n));

end

% Flip the backward trajectory
bwd_traj=fliplr(bwd_traj);

% Compute Re(<a|b>) fidelity
fidelity=real(rho_targ'*fwd_traj(:,end));

% Compute gradient
if nargout>1

    % Preallocate results
    grad=zeros(size(dt_grid));

    % Loop over control sequence
    parfor n=1:nsteps

        % Make evolution generator
        L=drift;
        for k=1:numel(controls)
            L=L+waveform(k,n)*controls{k};
        end

        % Compute fidelity derivative
        grad(n)=real(-1i*bwd_traj(:,n)'*L*fwd_traj(:,n));

    end

    % Convert time units
    grad=grad*time_unit;

end

end

% Consistency enforcement
function grumble(spin_system,drift,controls,waveform,dt_grid,rho_init,rho_targ)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('optimal control module requires Lioville space formalism.');
end
if (~isnumeric(rho_init))||(~iscolumn(rho_init))
    error('rho_init must be a column vector.');
end
if (~isnumeric(rho_targ))||(~iscolumn(rho_targ))
    error('rho_targ must be a column vector.');
end
if (~isnumeric(drift))||(size(drift,1)~=size(drift,2))
    error('drift must be a square matrix.');
end
if (size(drift,1)~=numel(rho_init))||(size(drift,1)~=numel(rho_targ))
    error('dimensions of drift, rho_init, and rho_targ must be consistent.');
end
if ~iscell(controls)
    error('controls must be a cell array of square matrices.');
end
for n=1:numel(controls)
    if (~isnumeric(controls{n}))||...
       (size(controls{n},1)~=size(controls{n},2))||...
       (size(controls{n},1)~=size(drift,1))
        error('control operators must have the same size as drift operators.')
    end
end
if (~isnumeric(waveform))||(~isreal(waveform))
    error('waveform must be a real numeric array.');
end
if size(waveform,1)~=numel(controls)
    error('number of waveform rows must be equal to the number of controls.');
end
if (~isnumeric(dt_grid))||(~isreal(dt_grid))||(~iscolumn(dt_grid))
    error('dt_grid must be a real row vector.');
end
if numel(dt_grid)~=size(waveform,2)
    error('numel(dt_grid) must be equal to the number waveform columns.');
end
end

% When I was twenty I became apprenticed to an old master cabinetmaker in
% Vienna [...] Once he told me that he had worked for many years on vari-
% ous models of a perpetual motion machine, adding musingly: "They say you
% can't make it; but once it's been made they'll talk different!"
%
% Karl Popper, "Unended Quest"

