% Emulates the effect of a gradient pulse on the sample average density
% matrix using Edwards formalism. It is assumed that the effect of dif-
% fusion is negligible, that the gradient is linear, and that it is an-
% tisymmetric about the middle of the sample. Syntax: 
%
%       rho=grad_pulse(spin_system,rho,g_amp,s_len,g_dur,s_fac)
%
% Parameters:
%   
%         rho   - spin system state vector
%
%           L   - system Liouvillian
%
%       g_amp   - gradient amplitude, Gauss/cm
%
%       s_len   - sample length, cm
%
%       g_dur   - gradient pulse duration, seconds
%
%       s_fac   - gradient shape factor, use 1 for
%                 square gradient pulses
%
% Outputs:
%
%      rho - spin system state vector, integrated over
%            the spatial coordinate
%
% Note: the function integrates over sample coordinates - subsequent gra-
%       dient pulses would not refocus the magnetization that it has de-
%       focused. To simulate a gradient sandwich, use grad_sandw.m func-
%       tion. More information on the subject is available in Luke's pa-
%       per (http://dx.doi.org/10.1016/j.jmr.2014.01.011).
%
% Note: this function is OK for standalone crusher gradients; for more
%       sophisticated gradient work, use the imaging context.
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=grad_pulse.m>

function rho=grad_pulse(spin_system,L,rho,g_amp,s_len,g_dur,s_fac)

% Check consistency
grumble(spin_system,L,rho,g_amp,s_len,g_dur,s_fac)

% Inform the user
report(spin_system,'computing the effect of a gradient pulse...');

% Get effective gradient operator (WARNING: shifts are ignored)
G=1e-4*s_fac*g_amp*s_len*g_dur*carrier(spin_system,'all')/spin_system.inter.magnet;

% Check approximation applicability
if cheap_norm(L*G-G*L)>1e-6
    error('the Liouvillian does not commute with the gradient operator.');
end

% Run background Liouvillian evolution
rho=evolution(spin_system,L,[],rho,g_dur,1,'final');

% Account for rescaling of sample integral into [0,1]
rho=evolution(spin_system,G,[],rho,-1/2,1,'final');

% Compute sample integral using auxiliary matrix method
aux_mat=[0*G, 1i*speye(size(G)); 0*G, G]; aux_rho=[0*rho; rho];
aux_rho=evolution(spin_system,aux_mat,[],aux_rho,1,1,'final');

% Map back into the state vector space
rho=aux_rho(1:end/2,:);

% Inform the user
report(spin_system,'gradient propagation done.')

end

% Consistency enforcement
function grumble(spin_system,L,rho,g_amp,s_len,g_dur,s_fac)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('this function is only applicable in Liouville space.');
end
if (~isnumeric(L))||(~isnumeric(rho))||(~isnumeric(g_amp))||...
   (~isnumeric(s_len))||(~isnumeric(g_dur))||(~isnumeric(s_fac))
   error('all arguments apart from spin_system must be numeric.');
end
if (~isreal(g_amp))||(~isscalar(g_amp))
    error('g_amp must be a real scalar.');
end
if (~isreal(s_len))||(~isscalar(s_len))||(s_len<=0)
    error('s_len must be a positive real number.');
end
if (~isreal(g_dur))||(~isscalar(g_dur))||(g_dur<0)
    error('g_dur must be a real non-negative scalar.');
end
if (~isreal(s_fac))||(~isscalar(s_fac))||(s_fac<0)
    error('s_fac must be a real non-negative scalar.');
end
end

% I cannot understand why we idle discussing religion. If we are honest - and scientists
% have to be - we must admit that religion is a jumble of false assertions, with no basis
% in reality. The very idea of God is a product of the human imagination. It is quite
% understandable why primitive people, who were so much more exposed to the overpowering
% forces of nature than we are today, should have personified these forces in fear and
% trembling. But nowadays, when we understand so many natural processes, we have no need
% for such solutions. I can't for the life of me see how the postulate of an Almighty God
% helps us in any way. What I do see is that this assumption leads to such unproductive
% questions as why God allows so much misery and injustice, the exploitation of the poor
% by the rich and all the other horrors He might have prevented. If religion is still
% being taught, it is by no means because its ideas still convince us, but simply because
% some of us want to keep the lower classes quiet. Quiet people are much easier to govern
% than clamorous and dissatisfied ones. They are also much easier to exploit. 
%
% Paul Dirac

