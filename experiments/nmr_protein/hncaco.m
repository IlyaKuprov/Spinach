% Protein-specific HN(CA)CO experiment (Figure 7.41 of "Protein NMR 
% Spectroscopy", 2nd edition) using pre-set values of J-couplings
% used in the magnetisation transfer stages. The simulation uses 
% the bidirectional propagation method described in
%
%           http://dx.doi.org/10.1016/j.jmr.2014.04.002
%
% The sequence is hard-wired to work on 1H,13C,15N proteins and uses
% PDB labels to select spins that will be affected by otherwise ideal
% pulses. F1 is 15N, F2 is 13C, F3 is 1H. Syntax:
%
%               fid=hnca(spin_system,parameters,H,R,K)
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
% m.walker@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=hncaco.m>

function fid=hncaco(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Coherent evolution timesteps
t1.nsteps=parameters.npoints(1); t1.timestep=1./parameters.sweep(1);
t2.nsteps=parameters.npoints(2); t2.timestep=1./parameters.sweep(2);
t3.nsteps=parameters.npoints(3); t3.timestep=1./parameters.sweep(3);

% J-coupling evolution time
tau=abs(1/(4*parameters.J_nh));
delta=abs(1/(2*parameters.J_nh));

% Pulse operators all protons
Hp=operator(spin_system,'L+','1H');
Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i;

% Spin indices and pulse operator - only on CA carbons
CAs=strcmp('CA',spin_system.comp.labels); 
CAp=operator(spin_system,'L+',find(CAs));
CAx=(CAp+CAp')/2; CAy=(CAp-CAp')/2i;

% Pulse operators all N
Np=operator(spin_system,'L+','15N');
Nx=(Np+Np')/2; 

% Spin indices and pulse operator - only on CO carbons
COs=strcmp('C',spin_system.comp.labels);
COp=operator(spin_system,'L+',find(COs));
COx=(COp+COp')/2;

% Detection state is protons
coil=state(spin_system,'L+','1H');

%% Run the first half forward

% Start in N+ to emulate the INEPT block
rho=state(spin_system,'L+','15N');

% Decouple protons
L_decH=decouple(spin_system,L,[],{'1H'});

% First half of t1 evolution
rho_stack=evolution(spin_system,L_decH,[],rho,t1.timestep/2,t1.nsteps-1,'trajectory');

% Inversion pulse on CO
rho_stack=step(spin_system,COx,rho_stack,pi);

% T/2 evolution
rho_stack=evolution(spin_system,L_decH,[],rho_stack,parameters.T/2,1,'final');

% Inversion pulses on N and CA
rho_stack=step(spin_system,Nx+CAx,rho_stack,pi);

% T/2 evolution
rho_stack=evolution(spin_system,L_decH,[],rho_stack,parameters.T/2,1,'final');

% -t1/2 evolution
rho_stack=evolution(spin_system,L_decH,[],rho_stack,-t1.timestep/2,t1.nsteps-1,'refocus');
                
% Pulses on 15N and 13CA
rho_stack=step(spin_system,Nx+CAx,rho_stack,pi/2);

% delta2 evolution
rho_stack=evolution(spin_system,L_decH,[],rho_stack,parameters.delta2,1,'final');

% Inversion pulses on CO and CA
rho_stack=step(spin_system,COx+CAx,rho_stack,pi);

% delta2 evolution
rho_stack=evolution(spin_system,L_decH,[],rho_stack,parameters.delta2,1,'final');

% Pulse on 13CO, y pulse on 13CA
rho_stack=step(spin_system,CAy+COx,rho_stack,pi/2);

%% Run the second half backward

% Get decoupled evolution generator
[L_decN,~]=decouple(spin_system,L,[],{'15N'});

% Detection on 1H
coil_stack=evolution(spin_system,L_decN,[],coil,...
                     -t3.timestep,t3.nsteps-1,'trajectory');

% tau evolution
coil_stack=evolution(spin_system,L,[],coil_stack,-tau,1,'final');

% Inversion pulses on 1H and 15N
coil_stack=step(spin_system,Hx+Nx,coil_stack,-pi);

% tau evolution
coil_stack=evolution(spin_system,L,[],coil_stack,-tau,1,'final');

% Pulses on 1H and 15N
coil_stack=step(spin_system,Hy+Nx,coil_stack,-pi/2);

% delta evolution
coil_stack=evolution(spin_system,L,[],coil_stack,-delta,1,'final');

% T/2 - delta evolution
coil_stack=evolution(spin_system,L_decH,[],coil_stack,-parameters.T/2+delta,1,'final');

% Inversion pulses on 13CA and 15N
coil_stack=step(spin_system,CAx+Nx,coil_stack,-pi);

% T/2 evolution
coil_stack=evolution(spin_system,L_decH,[],coil_stack,-parameters.T/2,1,'final');

% Pulse on 15N, y pulse on 13CA
coil_stack=step(spin_system,CAy+Nx,coil_stack,-pi/2);

% delta2 evolution
coil_stack=evolution(spin_system,L_decH,[],coil_stack,-parameters.delta2,1,'final');

% Inversion pulses on 13CA and 13CO
coil_stack=step(spin_system,CAx+COx,coil_stack,-pi);

% delta2 evolution
coil_stack=evolution(spin_system,L_decH,[],coil_stack,-parameters.delta2,1,'final');

% Pulses on 13CA and 13CO
coil_stack=step(spin_system,CAx+COx,coil_stack,-pi/2);

%% Stitch the halves

% Coherence selection
rho_stack=coherence(spin_system,rho_stack,{{'13C',+1}});
coil_stack=coherence(spin_system,coil_stack,{{'13C',+1}});

% Stitching
fid=stitch(spin_system,L_decH,rho_stack,coil_stack,{Nx+CAx},{pi},t1,t2,t3);

% Dimension reordering
fid=permute(fid,[3 2 1]);

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalisms.');
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
if ~isfield(parameters,'J_nh')
    error('scalar coupling should be specified in parameters.J_nh variable.');
elseif numel(parameters.J_nh)~=1
    error('parameters.J_nh array should have exactly one element.');
end
end

% Censorious puritan colleague: "How did you two meet?"
% IK's Tinder date: "In a church, over a game of scrabble."

