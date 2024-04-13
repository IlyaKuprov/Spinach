% Protein-specific H(CA)NH experiment (Figure 7.37 of "Protein NMR 
% Spectroscopy", 2nd edition) using pre-set values of J-couplings
% used in the magnetisation transfer stages. The simulation uses 
% the bidirectional propagation method described in
%
%              http://dx.doi.org/10.1016/j.jmr.2014.04.002
%
% The sequence is hard-wired to work on 1H,13C,15N proteins and uses
% PDB labels to select spins that will be affected by otherwise ideal
% pulses. F1 is 1H, F2 is 15N, F3 is 1H. Syntax:
%
%            fid=hcanh(spin_system,parameters,H,R,K)
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
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=hcanh.m>

function fid=hcanh(spin_system,parameters,H,R,K)

% Consistency check
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Coherent evolution timesteps
t1.nsteps=parameters.npoints(1); t1.timestep=1./parameters.sweep(1);
t2.nsteps=parameters.npoints(2); t2.timestep=1./parameters.sweep(2);
t3.nsteps=parameters.npoints(3); t3.timestep=1./parameters.sweep(3);

% Hard-coded J-couplings (Hz) and delays (seconds)
J_ch=140; J_nh=92;        % Textbook
tau1=abs(1/(4*J_ch));     % Sequence description
tau2=abs(1/(4*J_nh));     % Sequence description
delta1=abs(1/(4*J_ch));   % Sequence description
delta2=12.5e-3;           % Sequence description
delta3=abs(1/(4*J_nh));   % Sequence description
delta4=23.0e-3;           % Sequence description

% Initial condition - CA protons
if ~isfield(parameters,'rho0')
    HAs=ismember(spin_system.comp.labels,{'HA','HA1','HA2','HA3'}); 
    parameters.rho0=state(spin_system,'Lz',find(HAs),'cheap');
end

% Detection state - NH protons
if ~isfield(parameters,'coil')
    HNs=ismember(spin_system.comp.labels,{'H'}); 
    parameters.coil=state(spin_system,'L+',find(HNs),'cheap');
end

% Pulse operators all protons
Hp=operator(spin_system,'L+',parameters.spins{1});
Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i; 

% Spin indices and pulse operator - only on CA carbons
CAs=strcmp('CA',spin_system.comp.labels); 
CAp=operator(spin_system,'L+',find(CAs));
CAx=(CAp+CAp')/2;

% Pulse operators all N
Np=operator(spin_system,'L+','15N');
Nx=(Np+Np')/2; 

% Spin indices on CO carbons
COs=strcmp('C',spin_system.comp.labels);

%% Forward sim from rho0 up to t2 period 

% Pulse on 1H 
rho=step(spin_system,Hx,parameters.rho0,pi/2);

% Coherence selection
rho=coherence(spin_system,rho,{{'1H',1}});

% tau1 evolution
rho=evolution(spin_system,L,[],rho,tau1,1,'final');

% t1 evolution
rho_stack=evolution(spin_system,L,[],rho,t1.timestep/2,t1.nsteps-1,'trajectory');

% Inversion pulse on 13CA
rho_stack=step(spin_system,CAx,rho_stack,pi);

% t1 rest of the evolution
rho_stack=evolution(spin_system,L,[],rho_stack,t1.timestep/2,t1.nsteps-1,'refocus');    

% Inversion pulse on 1H
rho_stack=step(spin_system,Hx,rho_stack,pi);

% tau1 evolution
rho_stack=evolution(spin_system,L,[],rho_stack,tau1,1,'final');

% Pulse on 13CA, y pulse on 1H
rho_stack=step(spin_system,Hy+CAx,rho_stack,pi/2);

% Decoupling of 13CO
[L_decCO,rho_stack]=decouple(spin_system,L,rho_stack,find(COs));

% delta1 evolution
rho_stack=evolution(spin_system,L_decCO,[],rho_stack,delta1,1,'final');

% Inversion pulse on 1H
rho_stack=step(spin_system,Hx,rho_stack,pi);

% delta2-delta1 evolution
rho_stack=evolution(spin_system,L_decCO,[],rho_stack,delta2-delta1,1,'final');

% Inversion pulses on 13CA and 15N
rho_stack=step(spin_system,CAx+Nx,rho_stack,pi);

% delta2 evolution
rho_stack=evolution(spin_system,L_decCO,[],rho_stack,delta2,1,'final');

% Pulses on 13CA and 15N
rho_stack=step(spin_system,CAx+Nx,rho_stack,pi/2);

% delta3 evolution
rho_stack=evolution(spin_system,L_decCO,[],rho_stack,delta3,1,'final');

%% Backward sim from coil up to t2 period 

% Get decoupled evolution generator
[L_decCA,~]=decouple(spin_system,L,[],find(CAs));

% Detection on 1H
coil_stack=evolution(spin_system,L_decCA,[],parameters.coil,...
                     -t3.timestep,t3.nsteps-1,'trajectory');
                 
% Select single quantum coherence
coil_stack=coherence(spin_system,coil_stack,{{parameters.spins{1},1}});                 

% tau2 evolution
coil_stack=evolution(spin_system,L,[],coil_stack,-tau2,1,'final');

% Inversion pulses on 1H and 15N
coil_stack=step(spin_system,Hx+Nx,coil_stack,-pi);

% tau2 evolution
coil_stack=evolution(spin_system,L,[],coil_stack,-tau2,1,'final');

% Pulses on 1H and 15N
coil_stack=step(spin_system,Hx+Nx,coil_stack,-pi/2);

% delta4 evolution
coil_stack=evolution(spin_system,L_decCO,[],coil_stack,-delta4,1,'final');

% Inversion pulse on 15N
coil_stack=step(spin_system,Nx,coil_stack,-pi);

%% Stitch the halves

% Coherence selection
rho_stack=coherence(spin_system,rho_stack,{{'15N',+1}});
coil_stack=coherence(spin_system,coil_stack,{{'15N',+1}});

% Stitching
fid=stitch(spin_system,L_decCO,rho_stack,coil_stack,{Hx,L_decCO,CAx},...
           {pi,abs(delta4-delta3),pi},t1,t2,t3);

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
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=3
    error('parameters.sweep array should have exactly three elements.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
end
if (~isnumeric(parameters.sweep))||(~isvector(parameters.sweep))||...
   (~isreal(parameters.sweep))||(numel(parameters.sweep)~=3)
    error('parameters.sweep must be a vector of three real numbers.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(~isvector(parameters.npoints))||...
   (~isreal(parameters.npoints))||(numel(parameters.npoints)~=3)||...
    any(parameters.npoints<1)||any(mod(parameters.npoints,1)~=0)
    error('parameters.npoints must be a vector of three positive integers.');
end
end

% The running craze is a symptom of our deplorable age [...] Jogging
% is not only undignified but absurd. It is a confession that people
% feel that they lead displeasingly unhealthy lives, but are not pre-
% pared to do anything preventative, rather than remedial, about it.
% The answer for someone who thinks that he is overweight is to eat 
% less for a while, not to leap around at unseemly exercises. And 
% the way to eat less is, simply, to eat less.
%
% Geoffrey Wheatcroft

