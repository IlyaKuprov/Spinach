% Double-quantum version of the 13C-detected 14N-13C MAS 2D correlation
% experiment described by Jarvis, Haies, Williamson and Carravetta in
%
%                 http://dx.doi.org/10.1039/c3cp50787d
%
% Parameters:
%
%        parameters.rf_pwr   - RF power on 14N, Hz
%
%        parameters.rf_dur   - RF pulse duration on 14N
%
% Outputs:
%
%        fid.sin
%        fid.cos             - sine and cosine components
%                              of the States quadrature
%
% i.kuprov@soton.ac.uk
% p.t.williamson@soton.ac.uk
% m.carravetta@soton.ac.uk
% j.jarvis@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=cn2d_dq.m>

function fid=cn2d_dq(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose the Liouvillian
L=H+1i*R+1i*K;

% Get evolution timesteps
timestep=1./parameters.sweep;

% Get pulse operators
Cp=operator(spin_system,'L+',parameters.spins{2});
Np=operator(spin_system,'L+',parameters.spins{1});
Cx=(Cp+Cp')/2;  Cx=kron(speye(parameters.spc_dim),Cx);
Nx=(Np+Np')/2;  Nx=kron(speye(parameters.spc_dim),Nx);
Ny=(Np-Np')/2i; Ny=kron(speye(parameters.spc_dim),Ny);

% Apply 14N pulse
rho_cos=evolution(spin_system,L+2*pi*parameters.rf_pwr*Nx,[],parameters.rho0,parameters.rf_dur,1,'final');
rho_sin=evolution(spin_system,L+2*pi*parameters.rf_pwr*(Nx+Ny)/sqrt(2),[],parameters.rho0,parameters.rf_dur,1,'final');

% Apply coherence selection
rho_cos=coherence(spin_system,rho_cos,{{parameters.spins{2},[-1 1]},{parameters.spins{1},[-2 2]}});
rho_sin=coherence(spin_system,rho_sin,{{parameters.spins{2},[-1 1]},{parameters.spins{1},[-2 2]}});

% Run F1 evolution and 13C pulse
rho_stack_cos=evolution(spin_system,L,[],rho_cos,timestep(1)/2,parameters.npoints(1)-1,'trajectory');
rho_stack_sin=evolution(spin_system,L,[],rho_sin,timestep(1)/2,parameters.npoints(1)-1,'trajectory');
rho_stack_cos=step(spin_system,Cx,rho_stack_cos,pi); rho_stack_sin=step(spin_system,Cx,rho_stack_sin,pi);
rho_stack_cos=evolution(spin_system,L,[],rho_stack_cos,timestep(1)/2,parameters.npoints(1)-1,'refocus');
rho_stack_sin=evolution(spin_system,L,[],rho_stack_sin,timestep(1)/2,parameters.npoints(1)-1,'refocus');

% Apply 14N pulse
rho_stack_cos=evolution(spin_system,L+2*pi*parameters.rf_pwr*Nx,[],rho_stack_cos,parameters.rf_dur,1,'final');
rho_stack_sin=evolution(spin_system,L+2*pi*parameters.rf_pwr*Nx,[],rho_stack_sin,parameters.rf_dur,1,'final');

% Run F2 evolution
fid.cos=evolution(spin_system,L,parameters.coil,rho_stack_cos,timestep(2),parameters.npoints(2)-1,'observable');
fid.sin=evolution(spin_system,L,parameters.coil,rho_stack_sin,timestep(2),parameters.npoints(2)-1,'observable');
    
end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalism.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=2
    error('parameters.sweep array should have exactly two elements.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=2
    error('parameters.spins cell array should have exactly two elements.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=2
    error('parameters.npoints array should have exactly two elements.');
end
end

% It takes a great deal of expertise to distinguish between the
% three "I"s: innovators, imitators, and idiots.
%
% Warren Buffett

