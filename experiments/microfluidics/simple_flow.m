% Simple forward evolution experiment for the microfluidics 
% module; trajectory is returned. Syntax:
%
%    traj=simple_flow(spin_system,parameters,H,R,K,~,F)
%
% This sequence must be called from the meshflow() context,
% which would provide H, R, K, G, and F. Because gradients
% are not being used, the G input is ignored.
%
% Parameters:
%
%    parameters.npoints   - number of points in 
%                           the trajectory
%
%    parameters.rho0      - initial state in Fokker-
%                           Planck space
%
%    parameters.dt        - trajectory time step
%
% Outputs:
%
%    traj - trajectory in the Fokker-Planck space
%
% Notes: to convert Fokker-Planck space trajectory into R3
%        or Liouville space, use fpl2phan and fpl2rho func-
%        tions.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=simple_flow.m>

function traj=simple_flow(spin_system,parameters,H,R,K,~,F)

% Check consistency
grumble(parameters,H,R,K,F);

% Compose Liouvillian
L=H+1i*F+1i*R+1i*K;

% Run the evolution and watch the coil state
traj=evolution(spin_system,L,[],parameters.rho0,parameters.dt,...
               parameters.npoints-1,'trajectory');

end

% Consistency enforcement
function grumble(parameters,H,R,K,F)
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||(~isnumeric(F))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))||(~ismatrix(F))
    error('H, R, K, F arguments must be matrices.');
end
if (~isnumeric(parameters.dt))||(~isscalar(parameters.dt))||...
   (~isreal(parameters.dt))||(~isfinite(parameters.dt))
    error('parameters.dt should be finite real scalar.');
end
if (~isnumeric(parameters.npoints))||(~isscalar(parameters.npoints))||...
   (~isreal(parameters.npoints))||(parameters.npoints<1)||...
   (mod(parameters.npoints,1)~=0)
    error('parameters.npoints should be a positive integer.');
end
end

% Everybody has a plan until they get 
% punched in the face.
%
% Mike Tyson

