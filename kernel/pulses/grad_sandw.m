% Emulates the effect of a gradient sandwich on the sample average density
% matrix using Edwards formalism. It is assumed that the effect of diffusi-
% on is negligible, that the gradients are linear, and that they are anti-
% symmetric about the middle of the sample. Syntax:
%
%     rho=grad_sandw(spin_system,L,rho,P,g_amps,s_len,g_durs,s_facs)
%
% Parameters:
%   
%              rho   - spin system state vector
%
%                L   - system Liouvillian
%
%                P   - total propagator for all events happening
%                      between the two gradients
%
%           g_amps   - row vector containing the amplitudes of
%                      the two gradients, Gauss/cm
%
%            s_len   - sample length, cm
%
%           g_durs   - row vector containing the durations of
%                      the two gradients, seconds
%
%           s_facs   - shape factors of the two gradients, use
%                      [1 1] for square gradient pulses
%
% Outputs:
%
%      rho - spin system state vector, integrated over
%            the spatial coordinate
%
% Note: the function integrates over sample coordinates - subsequent gra-
%       dient pulses would not refocus the magnetization that it has left
%       defocused. More information on the subject is available in Luke's
%       paper (http://dx.doi.org/10.1016/j.jmr.2014.01.011).
%
% Note: this function is OK for standalone gradient pairs; for more
%       sophisticated gradient work, use the imaging context.
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=grad_sandw.m>

function rho=grad_sandw(spin_system,L,rho,P,g_amps,s_len,g_durs,s_facs)

% Check consistency
grumble(spin_system,L,rho,P,g_amps,s_len,g_durs,s_facs);

% Inform the user
report(spin_system,'computing the effect of a gradient sandwich...');

% Get effective gradient operators (WARNING: shifts are ignored)
R=carrier(spin_system,'all')/spin_system.inter.magnet;
G1=1e-4*s_facs(1)*g_amps(1)*s_len*g_durs(1)*R;
G2=1e-4*s_facs(2)*g_amps(2)*s_len*g_durs(2)*R;

% Check approximation applicability
if (cheap_norm(L*G1-G1*L)>1e-6)||...
   (cheap_norm(L*G2-G2*L)>1e-6)
    error('the Liouvillian does not commute with gradient operators.');
end

% Run background Liouvillian evolution during first gradient
rho=evolution(spin_system,L,[],rho,g_durs(1),1,'final');

% Account for rescaling of sample integral into [0,1]
rho=evolution(spin_system,G1,[],rho,-1/2,1,'final');

% Evolution under the gradient pulses and P
aux_mat=[-G2, 1i*P; 0*G1, G1]; aux_rho=[0*rho; rho];
aux_rho=evolution(spin_system,aux_mat,[],aux_rho,1,1,'final');

% Map back into the state vector space
rho=aux_rho(1:end/2,:);

% Undo the evolution overshoot
rho=evolution(spin_system,G2,[],rho,1/2,1,'final');

% Run background Liouvillian evolution during second gradient
rho=evolution(spin_system,L,[],rho,g_durs(2),1,'final');

end

% Consistency enforcement
function grumble(spin_system,L,rho,P,g_amps,s_len,g_durs,s_facs)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('this function is only applicable in Liouville space.');
end
if (~isnumeric(L))||(~isnumeric(rho))||(~isnumeric(P))||(~isnumeric(g_amps))||...
   (~isnumeric(s_len))||(~isnumeric(g_durs))||(~isnumeric(s_facs))
   error('all arguments apart from spin_system must be numeric.');
end
if (~isreal(g_amps))||(numel(g_amps)~=2)
    error('g_amps vector must have two real elements.');
end
if (~isreal(s_len))||(~isscalar(s_len))||(s_len<=0)
    error('s_len must be a positive real number.');
end
if (~isreal(g_durs))||(numel(g_durs)~=2)||(any(g_durs<0))
    error('g_durs vector must have two real non-negative elements.');
end
if (~isreal(s_facs))||(numel(s_facs)~=2)||(any(s_facs<0))
    error('s_facs vector must have two real non-negative elements.');
end
end

% Morality is simply the attitude we adopt towards 
% people we personally dislike.
%
% Oscar Wilde

