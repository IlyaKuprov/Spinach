% Protein-specific HNCA experiment (Figure 7.31a of "Protein NMR 
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
%    fid - a structure with four fields: fid.pos_pos, fid.pos_neg,
%          fid.neg_pos, fid.neg_neg that are used in the subsequ-
%          ent States quadrature processing
%
% Note: spin labels must be set to PDB atom IDs ('CA', 'HA', etc.) in
%       sys.labels for this sequence to work properly.
%
% m.walker@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=hnca.m>

function fid=hnca(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Coherent evolution timesteps
t1.nsteps=parameters.npoints(1); t1.timestep=1./parameters.sweep(1);
t2.nsteps=parameters.npoints(2); t2.timestep=1./parameters.sweep(2);
t3.nsteps=parameters.npoints(3); t3.timestep=1./parameters.sweep(3);

% Hard-coded J-couplings and delays
J_nh=92; tau=abs(1/(4*J_nh));
J_nca=11.5; delta=abs(1/(4*J_nca));

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

% Pulse operators all protons
Hp=operator(spin_system,'L+','1H');
Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i; 

% Pulse operator on CA carbons
CAs=strcmp('CA',spin_system.comp.labels); 
CAp=operator(spin_system,'L+',find(CAs));
CAx=(CAp+CAp')/2;

% Pulse operators all nitrogens
Np=operator(spin_system,'L+','15N');
Nx=(Np+Np')/2; 

% Pulse operator on CO carbons
COs=strcmp('C',spin_system.comp.labels);
COp=operator(spin_system,'L+',find(COs));
COx=(COp+COp')/2;

%% Run the first half forward

% Pulse on 1H 
rho=step(spin_system,Hx,parameters.rho0,pi/2);

% Coherence selection
rho=coherence(spin_system,rho,{{'1H',1}});

% tau evolution
rho=evolution(spin_system,L,[],rho,tau,1,'final');

% Inversion pulses on 1H and 15N
rho=step(spin_system,Hx+Nx,rho,pi);

% tau evolution
rho=evolution(spin_system,L,[],rho,tau,1,'final');

% Pulse on 15N, y pulse on 1H
rho=step(spin_system,Hy+Nx,rho,pi/2);

% Coherence selection for States quadrature in F1
rho_pos=coherence(spin_system,rho,{{'15N',+1}});
rho_neg=coherence(spin_system,rho,{{'15N',-1}});

% t1 evolution
rho_stack_pos=evolution(spin_system,L,[],rho_pos,t1.timestep/2,t1.nsteps-1,'trajectory');
rho_stack_neg=evolution(spin_system,L,[],rho_neg,t1.timestep/2,t1.nsteps-1,'trajectory');

% Inversion pulses on 1H, 13CA and 13CO
rho_stack_pos=step(spin_system,Hx+CAx+COx,rho_stack_pos,pi);
rho_stack_neg=step(spin_system,Hx+CAx+COx,rho_stack_neg,pi);

% t1 rest of the evolution
rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,t1.timestep/2,t1.nsteps-1,'refocus');
rho_stack_neg=evolution(spin_system,L,[],rho_stack_neg,t1.timestep/2,t1.nsteps-1,'refocus');

% delta evolution
rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,delta,1,'final');
rho_stack_neg=evolution(spin_system,L,[],rho_stack_neg,delta,1,'final');

% Pulses on 1H and 13CA
rho_stack_pos=step(spin_system,Hx+CAx,rho_stack_pos,pi/2);
rho_stack_neg=step(spin_system,Hx+CAx,rho_stack_neg,pi/2);

%% Run the second half backward

% Get decoupled evolution generator
[L_dec,~]=decouple(spin_system,L,[],{'15N'});

% Detection on 1H
coil_stack=evolution(spin_system,L_dec',[],parameters.coil,...
                     -t3.timestep,t3.nsteps-1,'trajectory');
                 
% Select single quantum coherence
coil_stack=coherence(spin_system,coil_stack,{{'1H',1}});                 

% tau evolution
coil_stack=evolution(spin_system,L',[],coil_stack,-tau,1,'final');

% Inversion pulses on 1H and 15N
coil_stack=step(spin_system,Hx+Nx,coil_stack,-pi);

% tau evolution
coil_stack=evolution(spin_system,L',[],coil_stack,-tau,1,'final');

% Pulses on 1H and 15N
coil_stack=step(spin_system,Hx+Nx,coil_stack,-pi/2);

% delta evolution
coil_stack=evolution(spin_system,L',[],coil_stack,-delta,1,'final');

% Pulses on 1H and 13CA
coil_stack=step(spin_system,Hx+CAx,coil_stack,-pi/2);

%% Stitch the halves

% Coherence selection for States quadrature in F2
coil_stack_pos=coherence(spin_system,coil_stack,{{find(CAs),+1}});
coil_stack_neg=coherence(spin_system,coil_stack,{{find(CAs),-1}});

% Stitching
fid.pos_pos=stitch(spin_system,L,rho_stack_pos,coil_stack_pos,{Hx+Nx+COx},{pi},t1,t2,t3);
fid.pos_neg=stitch(spin_system,L,rho_stack_pos,coil_stack_neg,{Hx+Nx+COx},{pi},t1,t2,t3);
fid.neg_pos=stitch(spin_system,L,rho_stack_neg,coil_stack_pos,{Hx+Nx+COx},{pi},t1,t2,t3);
fid.neg_neg=stitch(spin_system,L,rho_stack_neg,coil_stack_neg,{Hx+Nx+COx},{pi},t1,t2,t3);

% Dimension reordering
fid_names=fieldnames(fid);
for name_idx=1:numel(fid_names)
    fid.(fid_names{name_idx})=permute(fid.(fid_names{name_idx}),[3 2 1]);
end

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
   (~isreal(parameters.sweep))||(numel(parameters.sweep)~=3)||...
   any(~isfinite(parameters.sweep))||any(parameters.sweep<=0)
    error('parameters.sweep must be a vector of three positive real numbers.');
end
end

% "I have lately made an Experiment in Electricity that I
%  desire never to repeat." 
%
% Benjamin Franklin, in 1750, about accidentally 
% electrocuting himself while trying to kill a
% turkey with electricity.

