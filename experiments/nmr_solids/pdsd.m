% A simplified model of the PDSD experiment.
%
% <https://spindynamics.org/wiki/index.php?title=pdsd.m>

function fid=pdsd(spin_system,parameters,H,R,K)

% Compose the Liouvillian and its decoupled version
L=H+1i*R+1i*K; LD=decouple(spin_system,L,[],{'1H'});

% Get evolution timesteps
timestep=1/parameters.sweep;

% Get pulse operators
Cx=operator(spin_system,'Lx','13C');
Cx=kron(speye(parameters.spc_dim),Cx);
Hx=operator(spin_system,'Lx','1H');
Hx=kron(speye(parameters.spc_dim),Hx);

% Skip CP and start with Lx on carbons
rho=state(spin_system,'Lx','13C','cheap');
rho=kron(ones(parameters.spc_dim,1),rho);

% Quadrature detection state
coil=state(spin_system,'L+','13C','cheap');
coil=kron(ones(parameters.spc_dim,1),coil);
coil=coil/parameters.spc_dim;

% F1 evolution with proton decoupling
rho_stack=krylov(spin_system,LD,[],rho,timestep,...
                 parameters.npoints(1)-1,'trajectory');

% Second pulse
rho_stack=step(spin_system,Cx,rho_stack,pi/2);

% Mixing time evolution with proton irradiation
rho_stack=step(spin_system,L+2*pi*parameters.rate*Hx,...
               rho_stack,parameters.tmix);

% Third pulse
rho_stack=step(spin_system,Cx,rho_stack,pi/2);

% F2 evolution and detection with proton decoupling
fid=krylov(spin_system,LD,coil,rho_stack,timestep,...
           parameters.npoints(2)-1,'observable');
   
end

% A language is a dialect with an army and a navy.
% 
% Max Weinreich

