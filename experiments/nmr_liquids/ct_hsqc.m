% Constant-time phase-sensitive HSQC pulse sequence. Syntax:
%
%              fid=ct_hsqc(spin_system,parameters,H,R,K)
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
% ilya.kuprov@weizmann.ac.il
% mrw1g16@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ct_hsqc.m>

function fid=ct_hsqc(spin_system,parameters,H,R,K)

% Consistency check
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Coherent evolution timesteps
timestep=1./parameters.sweep;

% J-coupling evolution time
tau=abs(1/(4*parameters.J));

% Initial condition
if ~isfield(parameters,'rho0')
    parameters.rho0=state(spin_system,'Lz',parameters.spins{2},'cheap');
end

% Detection state
if ~isfield(parameters,'coil')
    parameters.coil=state(spin_system,'L+',parameters.spins{2},'cheap');
end

% Pulse operators
Cp=operator(spin_system,'L+',parameters.spins{1});
Cx=(Cp+Cp')/2;
Hp=operator(spin_system,'L+',parameters.spins{2});
Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i;

% Pulse on I spin (1H)
rho=step(spin_system,Hx,parameters.rho0,pi/2);

% tau evolution
rho=evolution(spin_system,L,[],rho,tau,1,'final');

% Inversion pulses
rho=step(spin_system,Hx+Cx,rho,pi);

% tau evolution
rho=evolution(spin_system,L,[],rho,tau,1,'final');

% Y pulse on I spin 
rho=step(spin_system,Hy,rho,pi/2);

% Pulse on S spin (13C) with coherence selection
rho=step(spin_system,Cx,rho,pi/2)-step(spin_system,Cx,rho,-pi/2);

% Preallocate rho stack
rho_stack=zeros([size(rho,1) parameters.npoints(1)],'like',1i);

% Get the time grid for the CT period
CT=parameters.npoints(1)/parameters.sweep(1);
t1_grid=linspace(0,CT,parameters.npoints(1));

% Loop over the value of t1_grid
for n=1:parameters.npoints(1)
    
    % Run the delay
    rho_current=evolution(spin_system,L,[],rho,(CT-t1_grid(n))/2,1,'final');
    
    % Run the pulse on the S spin
    rho_current=step(spin_system,Cx,rho_current,pi);
    
    % Run the delay
    rho_current=evolution(spin_system,L,[],rho_current,CT/2,1,'final');
    
    % Run the pulse on the I spin
    rho_current=step(spin_system,Hx,rho_current,pi);
    
    % Run the delay
    rho_current=evolution(spin_system,L,[],rho_current,t1_grid(n)/2,1,'final');
    
    % Assign the stack element
    rho_stack(:,n)=rho_current;
    
end

% Coherence selection - SQC on 13C is transferred to 1H in remaining steps
rho_stack_pos=coherence(spin_system,rho_stack,{{parameters.spins{2}, 0},...
                                               {parameters.spins{1},+1}});
rho_stack_neg=coherence(spin_system,rho_stack,{{parameters.spins{2}, 0},...
                                               {parameters.spins{1},-1}});

% Pulses on both spins
rho_stack_pos=step(spin_system,Cx+Hx,rho_stack_pos,pi/2);
rho_stack_neg=step(spin_system,Cx+Hx,rho_stack_neg,pi/2);

% tau evolution
rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,tau,1,'final');
rho_stack_neg=evolution(spin_system,L,[],rho_stack_neg,tau,1,'final');

% Inversion pulses on both spins
rho_stack_pos=step(spin_system,Cx+Hx,rho_stack_pos,pi);
rho_stack_neg=step(spin_system,Cx+Hx,rho_stack_neg,pi);

% Delta evolution
rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,tau,1,'final');
rho_stack_neg=evolution(spin_system,L,[],rho_stack_neg,tau,1,'final');

% Decoupling of S spin
[L,rho_stack_pos]=decouple(spin_system,L,rho_stack_pos,parameters.decouple_f2);
[L,rho_stack_neg]=decouple(spin_system,L,rho_stack_neg,parameters.decouple_f2);

% Detection on I spin
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

% She spent an astonishing amount of time in attending lec-
% tures and demonstrations, distributing literature for the
% Junior Anti-Sex League, preparing banners for Hate Week,
% making collections for the savings campaign, and such-like
% activities. It paid, she said; it was camouflage. If you
% kept the small rules you could break the big ones.
%
% George Orwell, "1984"

