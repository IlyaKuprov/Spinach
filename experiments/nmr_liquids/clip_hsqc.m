% CLIP-HSQC pulse sequence. Syntax:
%
%             fid=clip_hsqc(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep              [F1 F2] sweep widths, Hz
%
%    parameters.npoints            [F1 F2] numbers of points
%
%    parameters.spins              {F1 F2} nuclei (e.g. '1H', '13C')
%
%    parameters.J                  active scalar coupling, Hz
%
%    H     - Hamiltonian matrix, received from context function
%
%    R     - relaxation superoperator, received from context function
%
%    K     - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid   - free induction decay. States quadrature detection is
%            used; the two components of the States signal are re-
%            turned in fid.pos and fid.neg
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=clip_hsqc.m>

function fid=clip_hsqc(spin_system,parameters,H,R,K)

% Consistency check
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Coherent evolution timesteps
timestep=1./parameters.sweep;

% J-coupling evolution time
delta=abs(1/(2*parameters.J));

% Initial state
rho=state(spin_system,'Lz',parameters.spins{2},'cheap');

% Detection state
coil=state(spin_system,'L+',parameters.spins{2},'cheap');

% Pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Lx_F1=(Lp+Lp')/2;
Lp=operator(spin_system,'L+',parameters.spins{2});
Lx_F2=(Lp+Lp')/2; Ly_F2=(Lp-Lp')/2i;

% First pulse
report(spin_system,'propagating initial state forwards in time...')
rho=step(spin_system,Lx_F2,rho,pi/2);

% Delta evolution
rho=evolution(spin_system,L,[],rho,delta/2,1,'final');

% Second pulse
rho=step(spin_system,Lx_F2+Lx_F1,rho,pi);

% Delta evolution
rho=evolution(spin_system,L,[],rho,delta/2,1,'final');

% Pulses
rho=step(spin_system,Ly_F2,rho,pi/2);
rho_px=step(spin_system,+Lx_F1,rho,pi/2);
rho_mx=step(spin_system,-Lx_F1,rho,pi/2);

% F1 evolution
rho_stack_px=evolution(spin_system,L,[],rho_px,timestep(1)/2,parameters.npoints(1)-1,'trajectory');
rho_stack_mx=evolution(spin_system,L,[],rho_mx,timestep(1)/2,parameters.npoints(1)-1,'trajectory');
rho_stack_px=step(spin_system,Lx_F2,rho_stack_px,pi);
rho_stack_mx=step(spin_system,Lx_F2,rho_stack_mx,pi);
rho_stack_px=evolution(spin_system,L,[],rho_stack_px,timestep(1)/2,parameters.npoints(1)-1,'refocus');
rho_stack_mx=evolution(spin_system,L,[],rho_stack_mx,timestep(1)/2,parameters.npoints(1)-1,'refocus');

% Coherence selection by the first gradient
rho_stack_px_P=coherence(spin_system,rho_stack_px,{{parameters.spins{2},0},{parameters.spins{1},+1}});
rho_stack_mx_P=coherence(spin_system,rho_stack_mx,{{parameters.spins{2},0},{parameters.spins{1},+1}});
rho_stack_px_N=coherence(spin_system,rho_stack_px,{{parameters.spins{2},0},{parameters.spins{1},-1}});
rho_stack_mx_N=coherence(spin_system,rho_stack_mx,{{parameters.spins{2},0},{parameters.spins{1},-1}});

% F2 evolution
report(spin_system,'propagating coil state backwards in time...')
coil_stack=evolution(spin_system,L',[],coil,-timestep(2),parameters.npoints(2)-1,'trajectory');

% Pulses
coil_stack_px=step(spin_system,+Lx_F1,coil_stack,-pi/2);
coil_stack_mx=step(spin_system,-Lx_F1,coil_stack,-pi/2);

% Delta evolution
coil_stack_px=evolution(spin_system,L',[],coil_stack_px,-delta/2,1,'final');
coil_stack_mx=evolution(spin_system,L',[],coil_stack_mx,-delta/2,1,'final');

% Coherence selection by the second gradient
coil_stack_px=coherence(spin_system,coil_stack_px,{{parameters.spins{1},0},{parameters.spins{2},+1}});
coil_stack_mx=coherence(spin_system,coil_stack_mx,{{parameters.spins{1},0},{parameters.spins{2},+1}});

% Pulses
coil_stack_px=step(spin_system,Ly_F2+Lx_F1,coil_stack_px,-pi);
coil_stack_mx=step(spin_system,Ly_F2+Lx_F1,coil_stack_mx,-pi);

% Delta evolution
coil_stack_px=evolution(spin_system,L',[],coil_stack_px,-delta/2,1,'final');
coil_stack_mx=evolution(spin_system,L',[],coil_stack_mx,-delta/2,1,'final');

% Phase cycled components combined
coil_stack_px=step(spin_system,-Lx_F1,coil_stack_px,-3*pi/2)...
             -step(spin_system,+Lx_F1,coil_stack_px,-3*pi/2);
coil_stack_mx=step(spin_system,+Lx_F1,coil_stack_mx,-3*pi/2)...
             -step(spin_system,-Lx_F1,coil_stack_mx,-3*pi/2);
% Pulses
coil_stack_px=step(spin_system,+Lx_F2,coil_stack_px,-pi/2);
coil_stack_mx=step(spin_system,+Lx_F2,coil_stack_mx,-pi/2);

% Detect by taking the inner product of the forward and backward calculations
report(spin_system,'combining forward and backward calculations...')
fid.pos=coil_stack_px'*rho_stack_px_P+coil_stack_mx'*rho_stack_mx_P;
fid.neg=coil_stack_px'*rho_stack_px_N+coil_stack_mx'*rho_stack_mx_N;

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
    error('clip_hsqc: sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=2
    error('clip_hsqc: parameters.sweep array should have exactly two elements.');
end
if ~isfield(parameters,'spins')
    error('clip_hsqc: working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=2
    error('clip_hsqc: parameters.spins cell array should have exactly two elements.');
end
if ~isfield(parameters,'npoints')
    error('clip_hsqc: number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=2
    error('clip_hsqc: parameters.npoints array should have exactly two elements.');
end
if ~isfield(parameters,'J')
    error('clip_hsqc: scalar coupling should be specified in parameters.J variable.');
elseif numel(parameters.J)~=1
    error('clip_hsqc: parameters.J array should have exactly one element.');
end
end

% Shortly after the end of the Second World War, Joseph Stalin
% received a secret report detailing the amorous adventures of
% one of his top wartime commanders. General Rokossovsky's wo-
% men, it was indignantly stated, numbered in the dozens. The
% report requested orders from the Supreme Commander on what 
% to do. Stalin wrote in big bold letters: "ENVY".

