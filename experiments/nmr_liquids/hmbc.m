% Magnitude-mode HMBC pulse sequence from 
%
%   A. Bax and M.F. Summers, J. Am. Chem. Soc., 108, 2093 (1986)
%
% Syntax:
%
%              fid=hmbc(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep          [F1 F2] sweep widths in each
%                              frequency direction, Hz
%
%    parameters.npoints        [F1 F2] numbers of points in
%                              each time direction
%
%    parameters.spins          {F1 F2} nuclei, e.g. {'15N','1H'}
%
%    parameters.J              primary scalar coupling, Hz
%
%    parameters.delta_b        delta_2 delay from the paper
%                              cited above; the authors recom-
%                              mend 60e-3 seconds
%
%    H - Hamiltonian matrix, received from context function
%
%    R - relaxation superoperator, received from context function
%
%    K - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid - free induction decay for magnitude mode processing
%
% Note: natural abundance experiments should make use of the iso-
%       tope dilution functionality. See dilute.m function.
%
% Bud Macaulay, Ilya Kuprov
%
% <https://spindynamics.org/wiki/index.php?title=hmqc.m>

function fid=hmbc(spin_system,parameters,H,R,K)

% Consistency check
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Get evolution time step
timestep=1./parameters.sweep;

% Fixed evolution times
delta_a=abs(1/(2*parameters.J));
delta_b=parameters.delta_b;

% Initial and detection states
rho0=state(spin_system,'Lz',parameters.spins{2},'cheap');
coil=state(spin_system,'L+',parameters.spins{2},'cheap');

% Operators
Cp=operator(spin_system,'L+',parameters.spins{1});
Hp=operator(spin_system,'L+',parameters.spins{2});
Cx=(Cp+Cp')/2; Hx=(Hp+Hp')/2; 

% First proton pulse
rho=step(spin_system,Hx,rho0,pi/2);

% First evolution period
rho=evolution(spin_system,L,[],rho,delta_a,1,'final');

% First carbon pulse
rho=step(spin_system,Cx,rho,+pi/2);

% Second evolution period
rho=evolution(spin_system,L,[],rho,delta_b,1,'final');

% Second carbon pulse
rho=step(spin_system,Cx,rho,+pi/2)-...
    step(spin_system,Cx,rho,-pi/2);

% % Coherence selection
rho=coherence(spin_system,rho,{{'13C',+1}});

% F1 evolution period, first half
rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2,...
                    parameters.npoints(1)-1,'trajectory');

% Proton decoupling pulse
rho_stack=step(spin_system,Hx,rho_stack,pi);

% F1 evolution period, second half
rho_stack=evolution(spin_system,L,[],rho_stack,timestep(1)/2,...
                    parameters.npoints(1)-1,'refocus');

% Third carbon pulse
rho_stack=step(spin_system,Cx,rho_stack,pi/2);
      
% Detection
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
if ~isfield(parameters,'J')
    error('scalar coupling should be specified in parameters.J variable.');
elseif numel(parameters.J)~=1
    error('parameters.J array should have exactly one element.');
end
if ~isfield(parameters,'delta_b')
    error('second delay should be specified in parameters.delta_b variable.');
elseif numel(parameters.delta_b)~=1
    error('parameters.delta_b array should have exactly one element.');
end
end

% There is a beast in man that needs to be 
% exercised, not exorcised.
%
% Anton Szandor LaVey

