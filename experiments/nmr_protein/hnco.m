% Phase-sensitive HNCO pulse sequence from 
%
%             http://dx.doi.org/10.1016/0022-2364(90)90333-5
%
% using the bidirectional propagation method described in
%
%              http://dx.doi.org/10.1016/j.jmr.2014.04.002
%
% The sequence is hard-wired to work on 1H,13C,15N proteins and uses
% PDB labels to select spins that will be affected by otherwise ideal
% pulses. F1 is N, F2 is CO, F3 is H. Syntax:
%
%               fid=hnco(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.npoints     - a vector of three integers giving the
%                             number of points in the three temporal
%                             dimensions, ordered as [t1 t2 t3].
%
%    parameters.sweep       - a vector of three real numbers giving
%                             the sweep widths in the three frequen-
%                             cy dimensions, ordered as [f1 f2 f3].
%
%    parameters.tau         - the three delays required for the ope-
%                             ration of the sequence (see the paper)
%                             in seconds. Reasonable values are
%                             [2.25e-3, 14e-3, 4e-3]
%
%    parameters.f1_decouple - logical switch controlling proton de-
%                             coupling during the T1 period.
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
% <https://spindynamics.org/wiki/index.php?title=hnco.m>

function fid=hnco(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Evolution time discretization
t1.nsteps=parameters.npoints(1); t1.timestep=1./parameters.sweep(1);
t2.nsteps=parameters.npoints(2); t2.timestep=1./parameters.sweep(2);
t3.nsteps=parameters.npoints(3); t3.timestep=1./parameters.sweep(3);

% Initial and detection states
rho=state(spin_system,'Lz','1H','cheap');
coil=state(spin_system,'L+','1H','cheap');

% Spin indices
COs=strcmp('C',spin_system.comp.labels); 
CAs=strcmp('CA',spin_system.comp.labels);
Ns=strcmp('N',spin_system.comp.labels);

% Pulse operators
COp=operator(spin_system,'L+',find(COs));
CAp=operator(spin_system,'L+',find(CAs));
Hp=operator(spin_system,'L+','1H'); 
Np=operator(spin_system,'L+','15N'); 
Nx=(Np+Np')/2; COx=(COp+COp')/2; 
CAx=(CAp+CAp')/2; Hx=(Hp+Hp')/2;
Hy=(Hp-Hp')/2i;

% Run the first half forward
rho=step(spin_system,Hx,rho,pi/2);
rho=evolution(spin_system,L,[],rho,parameters.tau(1),1,'final');
rho=step(spin_system,Hx+Nx,rho,pi);
rho=evolution(spin_system,L,[],rho,parameters.tau(1),1,'final');
rho=step(spin_system,Hy+Nx,rho,pi/2);
rho=coherence(spin_system,rho,{{find(Ns),1}});
rho_stack=evolution(spin_system,L,[],rho,t1.timestep/2,t1.nsteps-1,'trajectory');
if parameters.f1_decouple
    rho_stack=step(spin_system,Hx+COx+CAx,rho_stack,pi);
else
    rho_stack=step(spin_system,COx+CAx,rho_stack,pi);
end
rho_stack=evolution(spin_system,L,[],rho_stack,t1.timestep/2,t1.nsteps-1,'refocus');
rho_stack=evolution(spin_system,L,[],rho_stack,parameters.tau(2),1,'final');
rho_stack=step(spin_system,CAx,rho_stack,pi);
rho_stack=evolution(spin_system,L,[],rho_stack,parameters.tau(3),1,'final');
rho_stack=step(spin_system,COx,rho_stack,pi/2);
rho_stack=coherence(spin_system,rho_stack,{{find(COs),1}});

% Run the second half backward
[L_dec,coil]=decouple(spin_system,L,coil,{'15N','13C'});
coil_stack=evolution(spin_system,L_dec',[],coil,-t3.timestep,t3.nsteps-1,'trajectory');
coil_stack=evolution(spin_system,L',[],coil_stack,-parameters.tau(1),1,'final');
coil_stack=step(spin_system,Hx+Nx,coil_stack,-pi);
coil_stack=evolution(spin_system,L',[],coil_stack,-parameters.tau(1),1,'final');
coil_stack=step(spin_system,Hx+Nx,coil_stack,-pi/2);
coil_stack=evolution(spin_system,L',[],coil_stack,-parameters.tau(2),1,'final');
coil_stack=step(spin_system,CAx,coil_stack,-pi);
coil_stack=evolution(spin_system,L',[],coil_stack,-parameters.tau(3),1,'final');
coil_stack=step(spin_system,COx,coil_stack,-pi/2);
coil_stack=coherence(spin_system,coil_stack,{{find(COs),1}});

% Stitch the halves
report(spin_system,'stitching forward and backward trajectories...');
fid=stitch(spin_system,L,rho_stack,coil_stack,{Nx+CAx},{pi},t1,t2,t3);

% Permute dimensions
fid=permute(fid,[3 2 1]);

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalisms.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'npoints')
    error('number of points must be specificed in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(~isvector(parameters.npoints))||(~isreal(parameters.npoints))||...
   (numel(parameters.npoints)~=3)||any(parameters.npoints<1)||any(mod(parameters.npoints,1)~=0)
    error('parameters.npoints must be a vector of three positive integers.');
end
if ~isfield(parameters,'sweep')
    error('sweep widths must be specificed in parameters.sweep variable.');
end
if (~isnumeric(parameters.sweep))||(~isvector(parameters.sweep))||(~isreal(parameters.sweep))||...
   (numel(parameters.sweep)~=3)
    error('parameters.sweep must be a vector of three real numbers.');
end
if ~isfield(parameters,'tau')
    error('sequence delays must be specificed in parameters.tau variable.');
end
if (~isnumeric(parameters.tau))||(~isvector(parameters.tau))||(~isreal(parameters.tau))||...
   (numel(parameters.tau)~=3)
    error('parameters.tau must be a vector of three real numbers.');
end
if ~isfield(parameters,'f1_decouple')
    error('F1 decoupling switch must be supplied in parameters.f1_decouple variable.');
end
if (~isnumeric(parameters.f1_decouple))||(~isscalar(parameters.f1_decouple))||...
   ((parameters.f1_decouple~=0)&&(parameters.f1_decouple~=1))
    error('parameters.J must be set to 1 or 0.');
end
end

% Envy is regarded by most people as a petty, superficial emotion
% and, therefore, it serves as a semihuman cover for so inhuman an
% emotion that those who feel it seldom dare admit it even to them-
% selves... That emotion is: hatred of the good for being the good.
%
% Ayn Rand

