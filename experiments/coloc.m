% COLOC NMR pulse sequence from 
%
%          https://doi.org/10.1016/0022-2364(84)90136-7
%
% implemented as shown in Figure 1b, without the dashed pulses during
% the delta(2) period. Delta(1) is calculated using the sweep width
% supplied by the user, delta(2) must be specified. Syntax:
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
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=coloc.m>

function fid=coloc(spin_system,parameters,H,R,K)

% Consistency check
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Sequence timing
timestep=1./parameters.sweep;

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
                     num2str(timestep(1)*(parameters.npoints(1)-1))]);

% First part of the echo block
rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2,...
                    parameters.npoints(1)-1,'trajectory');

% Echo pulse pair
rho_stack=step(spin_system,Hx,rho_stack,pi);
rho_stack=step(spin_system,Cx,rho_stack,pi);

% Second part of the echo block
rho_stack=fliplr(rho_stack);
rho_stack=evolution(spin_system,L,[],rho_stack,timestep(1)/2,...
                    parameters.npoints(1)-1,'refocus');
rho_stack=fliplr(rho_stack);

% Proton coherence selection
rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},-1}});

% 90 degree pulses on H and C
rho_stack=step(spin_system,Hx,rho_stack,pi/2);
rho_stack=step(spin_system,Cy,rho_stack,pi/2);

% Delta evolution
rho_stack=step(spin_system,L,rho_stack,parameters.delta2);

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
if ~isfield(parameters,'delta2')
    error('Delta2 delay should be specified in parameters.delta2 variable.');
elseif numel(parameters.delta2)~=1
    error('parameters.delta2 array should have exactly one element.');
end
end

% "NMR can be... annoyingly sensitive."
%
% Bernhard Bluemuch, 
% about water signals

