% A simplified model of the PDSD experiment.
%
% <https://spindynamics.org/wiki/index.php?title=pdsd.m>

function fid=pdsd_amp(spin_system,parameters,H,R,K)

% Compose the Liouvillian and its decoupled version
L=H+1i*R+1i*K; LD=decouple(spin_system,L,[],{'1H'});

% Get evolution timesteps
timestep=1/parameters.sweep;

% Get pulse operators
Cx=operator(spin_system,'Lx','13C');
Cx=kron(speye(parameters.spc_dim),Cx);
Hx=operator(spin_system,'Lx','1H');
Hx=kron(speye(parameters.spc_dim),Hx);

% Skip CP and start with Lz on carbons
rho=state(spin_system,'Lz','13C','cheap');
rho=kron(ones(parameters.spc_dim,1),rho);

% Quadrature detection state
coil=state(spin_system,'L+','13C','cheap');
coil=kron(ones(parameters.spc_dim,1),coil);
coil=coil/parameters.spc_dim;

% First pulse
rho=step(spin_system,Cx,rho,pi/2);

% F1 evolution
rho_stack=evolution(spin_system,L,[],rho,timestep,...
                    parameters.npoints(1)-1,'trajectory');
% Second pulse
rho_stack=step(spin_system,Cx,rho_stack,pi/2);

% Mixing time evolution with proton irradiation
rho_stack=evolution(spin_system,L+2*pi*parameters.rate*Hx,...
                    [],rho_stack,parameters.tmix,1,'final');

% Third pulse
rho_stack=step(spin_system,Cx,rho_stack,pi/2);

% Wipe the proton subspace (decoupling)
[~,rho_stack]=decouple(spin_system,[],rho_stack,{'1H'});
    
% F2 evolution and detection
fid=evolution(spin_system,LD,coil,rho_stack,timestep,...
              parameters.npoints(2)-1,'observable');
 
end

% A language is a dialect with an army and a navy.
% 
% Max Weinreich

