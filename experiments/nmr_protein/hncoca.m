% Magnitude-mode HN(CO)CA pulse sequence from 
%
%                http://dx.doi.org/10.1007/BF01874573
%
% using the bidirectional propagation method described in
%
%             http://dx.doi.org/10.1016/j.jmr.2014.04.002
%
% The sequence is hard-wired to work on 1H,13C,15N proteins and uses
% PDB labels to select spins that will be affected by otherwise ideal
% pulses. F1 is N, F2 is CA, F3 is H. Syntax:
%
%              fid=hncoca(spin_system,parameters,H,R,K)
%
% Parameters:
%
%      parameters.npoints   - a vector of three integers giving the
%                             number of points in the three temporal
%                             dimensions, ordered as [t1 t2 t3].
%
%      parameters.sweep     - a vector of three real numbers giving
%                             the sweep widths in the three frequen-
%                             cy dimensions, ordered as [f1 f2 f3].
%
%      parameters.tau       - the four delays required for the ope-
%                             ration of the sequence (see the paper)
%                             in seconds. Reasonable values are
%                             [2.25e-3, 2.75e-3, 8.00e-3, 7.00e-3]
%
%    H   - Hamiltonian matrix, received from context function
%
%    R   - relaxation superoperator, received from context function
%
%    K   - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid - three-dimensional free induction decay
%
% Note: spin labels must be set to PDB atom IDs ('CA', 'HA', etc.) in
%       sys.labels for this sequence to work properly.
% 
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=hncoca.m>

function fid=hncoca(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K)

% Compose Liouvillian
L=H+1i*R+1i*K;

% Evolution time discretization
t1.nsteps=parameters.npoints(1); t1.timestep=1/parameters.sweep(1);
t2.nsteps=parameters.npoints(2); t2.timestep=1/parameters.sweep(2);
t3.nsteps=parameters.npoints(3); t3.timestep=1/parameters.sweep(3);

% Hard-wired timings from the paper
tau1=2.25e-3; tau2=2.75e-3;
tau3=8.00e-3; tau4=7.00e-3;

% Initial condition - NH protons
if ~isfield(parameters,'rho0')
    HNs=ismember(spin_system.comp.labels,{'H'}); 
    parameters.rho0=state(spin_system,'Lz',find(HNs),'cheap');
end

% Detection state - NH protons
if ~isfield(parameters,'coil')
    HNs=ismember(spin_system.comp.labels,{'H'}); 
    parameters.coil=state(spin_system,'L+',find(HNs),'cheap');
end

% Spin indices
COs=strcmp('C',spin_system.comp.labels); 
CAs=strcmp('CA',spin_system.comp.labels);
Ns=strcmp('N',spin_system.comp.labels);

% Pulse operators
COp=operator(spin_system,'L+',find(COs));
CAp=operator(spin_system,'L+',find(CAs));
Hp=operator(spin_system,'L+','1H'); 
Np=operator(spin_system,'L+','15N');
Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i;
Nx=(Np+Np')/2; Ny=(Np-Np')/2i;
COx=(COp+COp')/2; CAx=(CAp+CAp')/2;

% Run the first half forward
rho=step(spin_system,Hx,parameters.rho0,pi/2);
rho=evolution(spin_system,L,[],rho,tau1,1,'final');
rho=step(spin_system,Hx+Ny,rho,pi);
rho=evolution(spin_system,L,[],rho,tau1,1,'final');
rho=step(spin_system,Hy+Nx,rho,pi/2);
rho=evolution(spin_system,L,[],rho,tau2,1,'final');
rho=coherence(spin_system,rho,{{find(Ns),1}});
rho_stack=evolution(spin_system,L,[],rho,t1.timestep/2,t1.nsteps-1,'trajectory');
rho_stack=step(spin_system,Hx,rho_stack,pi);
rho_stack=evolution(spin_system,L,[],rho_stack,t1.timestep/2,t1.nsteps-1,'refocus');
rho_stack=evolution(spin_system,L,[],rho_stack,tau3,1,'final');
rho_stack=step(spin_system,Nx+COx,rho_stack,pi);
rho_stack=evolution(spin_system,L,[],rho_stack,tau2,1,'final');
rho_stack=evolution(spin_system,L,[],rho_stack,tau3,1,'final');
rho_stack=step(spin_system,Nx+COx,rho_stack,pi/2);
rho_stack=evolution(spin_system,L,[],rho_stack,tau4,1,'final');
rho_stack=step(spin_system,CAx,rho_stack,pi/2);
rho_stack=coherence(spin_system,rho_stack,{{find(CAs),1}});

% Run the second half backward
[L_dec,coil]=decouple(spin_system,L,parameters.coil,find(strcmp('15N',spin_system.comp.isotopes)));
coil_stack=evolution(spin_system,L_dec',[],coil,-t3.timestep,t3.nsteps-1,'trajectory');
coil_stack=evolution(spin_system,L',[],coil_stack,-tau1,1,'final');
coil_stack=step(spin_system,Hx+Nx,coil_stack,-pi);
coil_stack=evolution(spin_system,L',[],coil_stack,-tau1,1,'final');
coil_stack=step(spin_system,Hx+Nx,coil_stack,-pi/2);
coil_stack=evolution(spin_system,L',[],coil_stack,-tau2,1,'final');
coil_stack=step(spin_system,Hx,coil_stack,-pi);
coil_stack=evolution(spin_system,L',[],coil_stack,-tau3,1,'final');
coil_stack=step(spin_system,Nx+COx,coil_stack,-pi);
coil_stack=evolution(spin_system,L',[],coil_stack,-tau2,1,'final');
coil_stack=evolution(spin_system,L',[],coil_stack,-tau3,1,'final');
coil_stack=step(spin_system,Nx+COx,coil_stack,-pi/2);
coil_stack=evolution(spin_system,L',[],coil_stack,-tau4,1,'final');
coil_stack=step(spin_system,CAx,coil_stack,-pi/2);
coil_stack=coherence(spin_system,coil_stack,{{find(CAs),1}});

% Stitch the halves
report(spin_system,'stitching forward and backward trajectories...');
fid=stitch(spin_system,L,rho_stack,coil_stack,{COx+Hx},{pi},t1,t2,t3);

% Permute dimensions
fid=permute(fid,[3 2 1]);

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalism.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'npoints')
    error('number of points must be specificed in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(~isvector(parameters.npoints))||...
   (~isreal(parameters.npoints))||(numel(parameters.npoints)~=3)||...
    any(parameters.npoints<1)||any(mod(parameters.npoints,1)~=0)
    error('parameters.npoints must be a vector of three positive integers.');
end
if ~isfield(parameters,'sweep')
    error('sweep widths must be specificed in parameters.sweep variable.');
end
if (~isnumeric(parameters.sweep))||(~isvector(parameters.sweep))||...
   (~isreal(parameters.sweep))||(numel(parameters.sweep)~=3)
    error('parameters.sweep must be a vector of three real numbers.');
end
end

% "This woman is headstrong, obstinate and dangerously self-opinionated."
%
% ICI personnel department assessment, 
% rejecting a job application from young 
% Margaret Thatcher in 1948.

