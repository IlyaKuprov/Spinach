% Phase-sensitive HSQC pulse sequence. Syntax:
%
%              fid=hsqc(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     parameters.sweep              [F1 F2] sweep widths, Hz
%
%     parameters.npoints            [F1 F2] numbers of points
%
%     parameters.spins              {F1 F2} nuclei (e.g. '13C','1H')
%
%     parameters.decouple_f2        nuclei to decouple in F2, e.g. 
%                                   {'15N','13C'}
%
%     parameters.decouple_f1        nuclei to decouple in F1, e.g. 
%                                   {'1H','13C'}
%
%     parameters.J                  working scalar coupling, Hz
%
%     H  - Hamiltonian matrix, received from context function
%
%     R  - relaxation superoperator, received from context function
%
%     K  - kinetics superoperator, received from context function
%
% Outputs:
%
%     fid.pos,fid.neg -  two components of the States quadrature 
%                        signal.
%
% Note: natural abundance simulations should make use of the isotope
%       dilution functionality. See dilute.m function.
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=hsqc.m>

function fid=hsqc(spin_system,parameters,H,R,K)

% Consistency check
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Coherent evolution timesteps
timestep=1./parameters.sweep;

% J-coupling evolution time
delta=abs(1/(2*parameters.J));

% Initial condition
if ~isfield(parameters,'rho0')
    parameters.rho0=state(spin_system,'Lz',parameters.spins{2},'cheap');
end

% Detection state
if ~isfield(parameters,'coil')
    parameters.coil=state(spin_system,'L+',parameters.spins{2},'cheap');
end

% Pulse operators
Lx_F1=operator(spin_system,'Lx',parameters.spins{1});
Lx_F2=operator(spin_system,'Lx',parameters.spins{2});
Ly_F2=operator(spin_system,'Ly',parameters.spins{2});

% Pulse on F2
rho=step(spin_system,Lx_F2,parameters.rho0,pi/2);

% Delta evolution
rho=evolution(spin_system,L,[],rho,delta/2,1,'final');

% Inversion pulses
rho=step(spin_system,Lx_F2+Lx_F1,rho,pi);

% Delta evolution
rho=evolution(spin_system,L,[],rho,delta/2,1,'final');

% Pulse on F2 
rho=step(spin_system,Ly_F2,rho,pi/2);

% Pulse on F1 with coherence selection
rho=step(spin_system,Lx_F1,rho,pi/2)-step(spin_system,Lx_F1,rho,-pi/2);

% F1 evolution
rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2,...
                    parameters.npoints(1)-1,'trajectory');
for n=1:numel(parameters.decouple_f1)
    Lp=operator(spin_system,'L+',parameters.decouple_f1{n});
    rho_stack=step(spin_system,(Lp+Lp')/2,rho_stack,pi);
end
rho_stack=evolution(spin_system,L,[],rho_stack,timestep(1)/2,...
                    parameters.npoints(1)-1,'refocus');

% Coherence selection
rho_stack_pos=coherence(spin_system,rho_stack,{{parameters.spins{2}, 0},...
                                               {parameters.spins{1},+1}});
rho_stack_neg=coherence(spin_system,rho_stack,{{parameters.spins{2}, 0},...
                                               {parameters.spins{1},-1}});

% Pulses on F1 and F2
rho_stack_pos=step(spin_system,Lx_F1+Lx_F2,rho_stack_pos,pi/2);
rho_stack_neg=step(spin_system,Lx_F1+Lx_F2,rho_stack_neg,pi/2);

% Delta evolution
rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,delta/2,1,'final');
rho_stack_neg=evolution(spin_system,L,[],rho_stack_neg,delta/2,1,'final');

% Pulses on F1 and F2
rho_stack_pos=step(spin_system,Lx_F1+Lx_F2,rho_stack_pos,pi);
rho_stack_neg=step(spin_system,Lx_F1+Lx_F2,rho_stack_neg,pi);

% Delta evolution
rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,delta/2,1,'final');
rho_stack_neg=evolution(spin_system,L,[],rho_stack_neg,delta/2,1,'final');

% Decoupling in F2
[L,rho_stack_pos]=decouple(spin_system,L,rho_stack_pos,parameters.decouple_f2);
[L,rho_stack_neg]=decouple(spin_system,L,rho_stack_neg,parameters.decouple_f2);

% Detect
fid.pos=evolution(spin_system,L,parameters.coil,rho_stack_pos,...
                  timestep(2),parameters.npoints(2)-1,'observable');
fid.neg=evolution(spin_system,L,parameters.coil,rho_stack_neg,...
                  timestep(2),parameters.npoints(2)-1,'observable');

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
if ~isfield(parameters,'J')
    error('scalar coupling should be specified in parameters.J variable.');
elseif numel(parameters.J)~=1
    error('parameters.J array should have exactly one element.');
end
end

% We do not hear the term "compassionate" applied to business executives or
% entrepreneurs, certainly not when they are engaged in their normal work.
% Yet in terms of results in the measurable form of jobs created, lives
% enriched, communities built, living standards raised, and poverty healed,
% a handful of capitalists has done infinitely more for mankind than all
% the self-serving politicians, academics, social workers, and religionists
% who march under the banner of "compassion".
%
% Nathaniel Branden

