% Multiple quantum correlation pulse sequence with refocusing. Syntax:
%
%            fid=mqs_refocus(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     parameters.sweep              [F1 F2] sweep widths (Hz)
%
%     parameters.npoints            [F1 F2] numbers of fid points
%
%     parameters.spins              {F1 F2} nuclei, in this case: 
%                                   {'1H','1H'} 
%
%     parameters.angle              flip angle for the final 
%                                   pulse, radians  
%
%     parameters.mqorder            coherence orders to select,
%                                   a two-element integer array
%
%     parameters.delay_1            first evolution delay, seconds
%
%     parameters.delay_2            second evolution delay, seconds
%
%     parameters.rho0               initial state
%
%     parameters.coil               detection state
%
%     H - Hamiltonian superoperator, provided by the context function
%
%     R - relaxation superoperator, provided by the context function
%
%     K - kinetics superoperator, provided by the context function
%
% Outputs:
%
%     fid - 2D free induction decay for amplitude-mode processing
% 
% mariagrazia.concilio@sjtu.edu.cn
%
% <https://spindynamics.org/wiki/index.php?title=mqs_refocus.m>

function fid=mqs_refocus(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Get pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Lx=(Lp+Lp')/2; Ly=(Lp-Lp')/2i;
    
% Apply the first 90 deg pulse
rho=step(spin_system,Lx,parameters.rho0,pi/2);
 
% Apply first evolution period
rho=evolution(spin_system,L,[],rho,parameters.delay_1,1,'final');

% Apply 180 deg pulse
rho=step(spin_system,Lx,rho,pi);

% Apply second evolution period
rho=evolution(spin_system,L,[],rho,parameters.delay_1,1,'final');

% Select operator
if mod(parameters.mqorder(1),2)==0
    
    % Apply the second 90 deg pulse about x
    rho=step(spin_system,Lx,rho,pi/2);    
    
else      
    % Apply the second 90 deg pulse about y
    rho=step(spin_system,Ly,rho,pi/2);        
end

% Coherence selection
rho=coherence(spin_system,rho,{{parameters.spins{1},parameters.mqorder(1)}});

% Run the F1 evolution
rho=evolution(spin_system,L,[],rho,1/parameters.sweep(1),...
              parameters.npoints(1)-1,'trajectory');

% Apply third pulse
rho=step(spin_system,Lx,rho,parameters.angle);

% Select coherence -1
rho=coherence(spin_system,rho,{{parameters.spins{1},parameters.mqorder(2)}});

% Apply first evolution period
rho=evolution(spin_system,L,[],rho,parameters.delay_2,1,'final');

% Apply refocussing pulse
rho=step(spin_system,Lx,rho,pi);

% Apply first evolution period
rho=evolution(spin_system,L,[],rho,parameters.delay_2,1,'final');

% Run the F2 evolution
fid=evolution(spin_system,L,parameters.coil,rho,1/parameters.sweep(2),...
                            parameters.npoints(2)-1,'observable');

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
    error('The sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=2
    error('The parameters.sweep array should have exactly two elements.');
elseif (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
       any(parameters.sweep<=0)
    error('parameters.sweep must contain two positive real numbers.');
end
if ~isfield(parameters,'spins')
    error('Two spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=2
    error('The parameters.spins cell array should have exactly two elements.');
elseif (~iscell(parameters.spins))||(~ischar(parameters.spins{1}))||...
       (~ischar(parameters.spins{2}))
    error('parameters.spins must be a two-element cell array of character strings.');
end
if ~isfield(parameters,'npoints')
    error('The number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=2
    error('The parameters.npoints array should have exactly two elements.');
elseif (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
       any(parameters.npoints<1)||any(mod(parameters.npoints,1)~=0)
    error('parameters.npoints must contain two positive integers.');
end
if ~isfield(parameters,'angle')
    error('pulse angle should be specified in parameters.angle variable.');
elseif numel(parameters.angle)~=1
    error('parameters.angle array should have exactly one element.');
elseif (~isnumeric(parameters.angle))||(~isreal(parameters.angle))
    error('parameters.angle must be a real scalar.');
end
if ~isfield(parameters,'mqorder')
    error('the multiple quantum coherence order should be specified in parameters.mqorder variable.');
elseif numel(parameters.mqorder)~=2
    error('parameters.mqorder array should have exactly two elements.');
elseif (~isnumeric(parameters.mqorder))||(~isreal(parameters.mqorder))||...
       any(mod(parameters.mqorder,1)~=0)
    error('parameters.mqorder must contain two integers.');
end
if ~isfield(parameters,'delay_1')
    error('first evolution delay should be specified in parameters.delay_1 variable.');
elseif numel(parameters.delay_1)~=1
    error('parameters.delay_1 array should have exactly one element.');
elseif (~isnumeric(parameters.delay_1))||(~isreal(parameters.delay_1))||...
       (parameters.delay_1<0)
    error('parameters.delay_1 must be a non-negative real scalar.');
end
if ~isfield(parameters,'delay_2')
    error('second evolution delay should be specified in parameters.delay_2 variable.');
elseif numel(parameters.delay_2)~=1
    error('parameters.delay_2 array should have exactly one element.');
elseif (~isnumeric(parameters.delay_2))||(~isreal(parameters.delay_2))||...
       (parameters.delay_2<0)
    error('parameters.delay_2 must be a non-negative real scalar.');
end
if ~isfield(parameters,'rho0')
    error('initial state should be specified in parameters.rho0 variable.');
elseif ~isnumeric(parameters.rho0)
    error('parameters.rho0 must be a numeric array.');
end
if ~isfield(parameters,'coil')
    error('detection state should be specified in parameters.coil variable.');
elseif ~isnumeric(parameters.coil)
    error('parameters.coil must be a numeric array.');
end
end

% Never underestimate the bandwidth of a station wagon full of
% tapes hurtling down the highway.
%
% Andrew Tanenbaum

