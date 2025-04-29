% A simplified model of the PDSD experiment.
%
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de
%
% <https://spindynamics.org/wiki/index.php?title=pdsd.m>

function fid=pdsd(spin_system,parameters,H,R,K)

% Compose the Liouvillian and its decoupled version
L=H+1i*R+1i*K; LD=decouple(spin_system,L,[],{'1H'});

% Get evolution timesteps
timestep=1/parameters.sweep;

% Get pulse operators
Cx=operator(spin_system,'Lx','13C');
Cy=operator(spin_system,'Ly','13C');
Cx=kron(speye(parameters.spc_dim),Cx);
Cy=kron(speye(parameters.spc_dim),Cy);
Hx=operator(spin_system,'Lx','1H');
Hx=kron(speye(parameters.spc_dim),Hx);

% Quadrature detection state
coil=state(spin_system,'L+','13C','cheap');
coil=kron(ones(parameters.spc_dim,1),coil);
coil=coil/parameters.spc_dim;

% Skip CP and start with Ly on carbons
rho=state(spin_system,'Ly','13C','cheap');
rho=kron(ones(parameters.spc_dim,1),rho);

% Phase cycle specification
Op2={Cx,Cy,Cx,Cy}; An2={+pi/2,+pi/2,-pi/2,-pi/2};
Op3={Cy,Cy,Cy,Cy}; An3={+pi/2,+pi/2,+pi/2,+pi/2};

% FID phase cycle
fids=cell(1,4);

% Phase cycle loop
for n=1:4

    % F1 evolution under proton decoupling
    rho_stack=evolution(spin_system,LD,[],rho,timestep,...
                        parameters.npoints(1)-1,'trajectory');
    % Second pulse
    rho_stack=step(spin_system,Op2{n},rho_stack,An2{n});

    % Mixing time evolution with proton irradiation
    rho_stack=evolution(spin_system,L+2*pi*parameters.rate*Hx,...
                        [],rho_stack,parameters.tmix,1,'final');

    % Third pulse
    rho_stack=step(spin_system,Op3{n},rho_stack,An3{n});

    % Wipe the proton subspace (decoupling)
    [~,rho_stack]=decouple(spin_system,[],rho_stack,{'1H'});
    
    % F2 evolution and detection under proton decoupling
    fids{n}=evolution(spin_system,LD,coil,rho_stack,timestep,...
                      parameters.npoints(2)-1,'observable');
                  
end

% Axial peak elimination
fid.cos=fids{1}-fids{3}; fid.sin=fids{2}-fids{4};
  
end

% A language is a dialect with an army and a navy.
% 
% Max Weinreich

