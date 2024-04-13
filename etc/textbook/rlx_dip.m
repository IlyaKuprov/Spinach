% Redfield theory expressions for dipolar relaxation and cross-
% relaxation rates, isotropic tumbling in liquid phase. Syntax:
%
%             [r1,r2,rx]=rlx_dip(B0,spins,dist,tau_c)
%
% Parameters:
%
%    B0     - magnet field, Tesla
%
%    spins  - the spins involved, e.g. {'1H','15N'}
%
%    dist   - inter-spin distance, Angstrom
%
%    tau_c  - rotational correlation time, seconds
%
% Outputs:
%
%    r1     - two longitudinal relaxation rates, Hz
%
%    r2     - two transverse relaxation rates, Hz
%
%    rx     - longitudinal cross-relaxation rate, Hz
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rlx_dip.m>

function [r1,r2,rx]=rlx_dip(B0,spins,dist,tau_c)

% Check consistency
grumble(B0,spins,dist,tau_c);

% Blicharsky invariant and rotational diffusion coefficient
[~,~,~,~,DD]=xyz2dd([0 0 0],[0 0 dist],spins{1},spins{2});
[~,DSqDD]=blinv(DD); r_dif_c=1/(6*tau_c);

% Multiplicities and spin-squares
[~,m_a]=spin(spins{1}); s_a=(m_a-1)/2; Ssq_a=s_a*(s_a+1);
[~,m_b]=spin(spins{2}); s_b=(m_b-1)/2; Ssq_b=s_b*(s_b+1);

% Zeeman frequencies
omega_a=spin(spins{1})*B0; omega_b=spin(spins{2})*B0;

% Textbook equation, R1 first spin
r1(1)=(2/27)*Ssq_b*DSqDD*(3*spden(2,r_dif_c,omega_a)+...
                          6*spden(2,r_dif_c,omega_a+omega_b)+...
                          1*spden(2,r_dif_c,omega_a-omega_b));
 
% Textbook equation, R1 second spin
r1(2)=(2/27)*Ssq_a*DSqDD*(3*spden(2,r_dif_c,omega_b)+...
                          6*spden(2,r_dif_c,omega_a+omega_b)+...
                          1*spden(2,r_dif_c,omega_a-omega_b));

% Textbook equation, R2 first spin
r2(1)=(1/27)*Ssq_b*DSqDD*(4*spden(2,r_dif_c,0)+...
                          3*spden(2,r_dif_c,omega_a)+...
                          6*spden(2,r_dif_c,omega_b)+...
                          6*spden(2,r_dif_c,omega_a+omega_b)+...
                          1*spden(2,r_dif_c,omega_a-omega_b));

% Textbook equation, R2 second spin
r2(2)=(1/27)*Ssq_a*DSqDD*(4*spden(2,r_dif_c,0)+...
                          3*spden(2,r_dif_c,omega_b)+...
                          6*spden(2,r_dif_c,omega_a)+...
                          6*spden(2,r_dif_c,omega_a+omega_b)+...
                          1*spden(2,r_dif_c,omega_a-omega_b));

% Textbook equation, longitudinal cross-relaxation rate
rx=(2/27)*sqrt(Ssq_a*Ssq_b)*DSqDD*(6*spden(2,r_dif_c,omega_a+omega_b)-...
                                   1*spden(2,r_dif_c,omega_a-omega_b));
    
end

% Consistency enforcement
function grumble(B0,spins,dist,tau_c)
if (~isnumeric(B0))||(~isreal(B0))||(~isscalar(B0))
    error('B0 must be a real number.');
end
if (~isnumeric(tau_c))||(~isreal(tau_c))||...
   (~isscalar(tau_c))||(tau_c<=0)
    error('tau_c must be a positive real number.');
end
if (~isnumeric(dist))||(~isreal(dist))||...
   (~isscalar(dist))||(dist<=0)
    error('dist must be a positive real number.');
end
if (~iscell(spins))||(numel(spins)~=2)||...
   (~ischar(spins{1}))||(~ischar(spins{2}))
    error('spins must be a cell array containing two character strings.');
end
end

% Provost Phelps [of Oriel College] was in the habit of taking a cold bath
% each morning. His butler would hear him saying, before he stepped into
% the icy waters: "Be a man, Phelps!" He was also the first Oxford head of
% house to entertain a puritanical new head of the nonconformist Mansfield
% College. When the port came round, this prig remarked: "I'd rather commit
% adultery". "Wouldn't we all?" - asked Phelps, before liberally refilling
% his own glass.
%
% A.N. Wilson

