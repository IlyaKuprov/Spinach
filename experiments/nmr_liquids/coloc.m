% COLOC NMR pulse sequence from 
%
%          https://doi.org/10.1016/0022-2364(84)90136-7
%
% implemented as shown in Figure 1b, without the dashed pulses during
% the delta(2) period. Delta(1) defaults to half of the maximum F1
% evolution time implied by the sweep width, delta(2) must be specified.
% Syntax:
%
%             fid=coloc(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     parameters.sweep         [F1 F2] sweep widths, Hz
%
%     parameters.npoints       [F1 F2] numbers of points
%
%     parameters.spins         {F1 F2} nuclei (e.g. '13C','1H')
%
%     parameters.delta2        COLOC delta2 (see the paper),
%                              typically 40e-3 seconds
%
%     parameters.delta1        optional COLOC delta1 delay, seconds;
%                              must be at least half of maximum t1
%
%     H  - Hamiltonian matrix, received from context function
%
%     R  - relaxation superoperator, received from context function
%
%     K  - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid - free induction decay for magnitude mode processing
%
% Note: natural abundance simulations should make use of the isotope
%       dilution functionality. See dilute.m function.
% 
% b.macaulay@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=coloc.m>

function fid=coloc(spin_system,parameters,H,R,K)

% Consistency check
grumble(spin_system,parameters,H,R,K);

% Set default delta1 delay
if ~isfield(parameters,'delta1')
    parameters.delta1=(parameters.npoints(1)-1)/(2*parameters.sweep(1));
end

% Compose Liouvillian
L=H+1i*R+1i*K;

% Sequence timing
timestep=1./parameters.sweep;
t1_grid=(0:(parameters.npoints(1)-1))*timestep(1);

% Initial state
rho=state(spin_system,'Lz',parameters.spins{1},'cheap');

% Detection state
coil=state(spin_system,'L+',parameters.spins{2},'cheap');

% Pulse operators
Hp=operator(spin_system,'L+',parameters.spins{1});
Cp=operator(spin_system,'L+',parameters.spins{2}); 
Hx=(Hp+Hp')/2; Cx=(Cp+Cp')/2; Cy=(Cp-Cp')/2i;

% 90 degree pulse on H
rho=step(spin_system,Hx,rho,pi/2);

% Report the delta(1)
report(spin_system,['COLOC delta(1), seconds: ' ...
                     num2str(parameters.delta1)]);

% Preallocate the evolution stack
rho_stack=zeros([size(rho,1) parameters.npoints(1)],'like',1i);

% Loop over the embedded-echo evolution time
for n=1:parameters.npoints(1)

    % First part of the echo block
    rho_current=evolution(spin_system,L,[],rho,...
                          t1_grid(n)/2,1,'final');

    % Echo pulse pair
    rho_current=step(spin_system,Hx,rho_current,pi);
    rho_current=step(spin_system,Cx,rho_current,pi);

    % Second part of the echo block
    rho_current=evolution(spin_system,L,[],rho_current,...
                          parameters.delta1-t1_grid(n)/2,1,'final');

    % Assign the stack element
    rho_stack(:,n)=rho_current;

end

% Proton coherence selection
rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},-1}});

% 90 degree pulses on H and C
rho_stack=step(spin_system,Hx,rho_stack,pi/2);
rho_stack=step(spin_system,Cy,rho_stack,pi/2);

% Delta evolution
rho_stack=evolution(spin_system,L,[],rho_stack,parameters.delta2,1,'final');

% Carbon coherence selection
rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{2},+1}});

% Proton decoupling
[L,rho_stack]=decouple(spin_system,L,rho_stack,parameters.spins(1));

% Detection on C
fid=evolution(spin_system,L,coil,rho_stack,timestep(2),...
              parameters.npoints(2)-1,'observable');

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
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=2
    error('parameters.sweep array should have exactly two elements.');
elseif (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
       any(~isfinite(parameters.sweep))||any(parameters.sweep<=0)
    error('parameters.sweep must contain two positive real numbers.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=2
    error('parameters.spins cell array should have exactly two elements.');
elseif (~iscell(parameters.spins))||(~ischar(parameters.spins{1}))||...
       (~ischar(parameters.spins{2}))
    error('parameters.spins must be a two-element cell array of character strings.');
elseif any(~ismember(parameters.spins,spin_system.comp.isotopes))
    error('parameters.spins contains isotopes that are not present in the system.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=2
    error('parameters.npoints array should have exactly two elements.');
elseif (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
       any(parameters.npoints<1)||any(mod(parameters.npoints,1)~=0)
    error('parameters.npoints must contain two positive integers.');
end
if ~isfield(parameters,'delta2')
    error('Delta2 delay should be specified in parameters.delta2 variable.');
elseif numel(parameters.delta2)~=1
    error('parameters.delta2 array should have exactly one element.');
elseif (~isnumeric(parameters.delta2))||(~isreal(parameters.delta2))||...
       (~isfinite(parameters.delta2))||(parameters.delta2<=0)
    error('parameters.delta2 must be a positive real scalar.');
end
if isfield(parameters,'delta1')
    max_t1=(parameters.npoints(1)-1)/parameters.sweep(1);
    if numel(parameters.delta1)~=1
        error('parameters.delta1 array should have exactly one element.');
    elseif (~isnumeric(parameters.delta1))||(~isreal(parameters.delta1))||...
           (~isfinite(parameters.delta1))||(parameters.delta1<0)
        error('parameters.delta1 must be a non-negative real scalar.');
    elseif parameters.delta1<max_t1/2
        error('parameters.delta1 must be at least half of the maximum F1 evolution time.');
    end
end
end

% "NMR can be... annoyingly sensitive."
%
% Bernhard Bluemuch, 
% about water signals

