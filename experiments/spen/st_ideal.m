% The ideal Stejskal-Tanner pulse sequence using the notation from
% Figure 1 in http://dx.doi.org/0.1002/cmr.a.21241 with no gaps be-
% tween pulse sequence events. Syntax:
%
%        inten=st_ideal(spin_system,parameters,H,R,K,G,F)
%
% This sequence must be called from the imaging() context, which
% would provide H, R, K, G, and F. Parameters:
%
%      parameters.spins      - working spin.
%
%      parameters.g_amp      - gradient amplitude, T/m
%
%      parameters.delta_sml  - the small delta parameter
%                              (see the figure)
%
%      parameters.delta_big  - the big delta parameter
%                              (see the figure)
%
% Outputs:
%
%      inten  - the absolute value of the first point in
%               the free induction decay; this number is 
%               proportional to the integral of the real
%               part of the correctly phased spectrum
%
%
% mariagrazia.concilio@sjtu.edu.cn
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=st_ideal.m>

function inten=st_ideal(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,F,G)

% Compose Liouvillian
L=H+F+1i*R+1i*K;

% Get pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i);

% Apply the excitation pulse
rho=step(spin_system,Ly,parameters.rho0,pi/2);

% Evolve under the first gradient
rho=evolution(spin_system,L+parameters.g_amp*G{1},[],rho,parameters.delta_sml,1,'final');

% Diffusion delay
rho=step(spin_system,L,rho,(parameters.delta_big-parameters.delta_sml)/2);

% Apply the refocusing pulse
rho=step(spin_system,Ly,rho,pi);

% Evolve to the second gradient
rho=step(spin_system,L,rho,(parameters.delta_big-parameters.delta_sml)/2);

% Evolve under the second gradient
rho=evolution(spin_system,L+parameters.g_amp*G{1},[],rho,parameters.delta_sml,1,'final');

% Intensity is the first point in the FID
inten=abs(parameters.coil'*rho);

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K,F,G)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalism.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~isnumeric(F))||(~ismatrix(H))||(~ismatrix(R))||...
   (~ismatrix(K))||(~ismatrix(F))
    error('H, R, K and F must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))||(~all(size(K)==size(F)))
    error('H, R, K and F matrices must have the same dimension.');
end
if ~iscell(G), error('G must be a cell array.'); end
if ~isfield(parameters,'g_amp')
    error('the gradient amplitude should be specified in parameters.g_amp variable.');
elseif numel(parameters.g_amp)~=1
    error('parameters.g_amp array should have exactly one element.');
end
if ~isfield(parameters,'delta_sml')
    error('the gradient duration should be specified in parameters.delta_sml variable.');
elseif numel(parameters.delta_sml)~=1
    error('parameters.delta_sml array should have exactly one element.');
end
if ~isfield(parameters,'delta_big')
    error('the diffusion delay should be specified in parameters.delta_sml variable.');
elseif numel(parameters.delta_sml)~=1
    error('parameters.delta_sml array should have exactly one element.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=1
    error('parameters.spins cell array should have exactly one element.');
end
end

% Whether it be movie stars, or athletes, or musicians, or CEOs of
% multibillion-dollar corporations who dropped out of school, popu-
% lar media often tells the story of the determined individual who
% pursues their dreams and beats the odds. [...] This creates a fal-
% se public perception that anyone can achieve great things if they
% have the ability and make the effort. The overwhelming majority
% of failures are not visible to the public eye, and only those who
% survive the selective pressures of their competitive environment
% are seen regularly.
%
% Wikipedia, about survivorship bias

