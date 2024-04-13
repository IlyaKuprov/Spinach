% Redfield theory expressions for hyperfine relaxation and cross-
% relaxation rates, isotropic tumbling in liquid phase. Syntax:
%
%             [r1,r2,rx]=rlx_hfc(B0,A,spins,tau_c)
%
% Parameters:
%
%    B0     - magnet field, Tesla
%
%    A      - 3x3 hyperfine coupling tensor,
%             not necessarily symmetric, rad/s
%
%    spins  - the spins involved, e.g. {'E','15N'},
%             one of those must be an electron
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
% <https://spindynamics.org/wiki/index.php?title=rlx_hfc.m>

function [r1,r2,rx]=rlx_hfc(B0,HFC,spins,tau_c)

% Check consistency
grumble(B0,HFC,spins,tau_c);

% Blicharsky invariants and diffusion coefficient
[LSqHFC,DSqHFC]=blinv(HFC); r_dif_c=1/(6*tau_c);

% Multiplicities and spin-squares
[~,m_a]=spin(spins{1}); s_a=(m_a-1)/2; Ssq_a=s_a*(s_a+1);
[~,m_b]=spin(spins{2}); s_b=(m_b-1)/2; Ssq_b=s_b*(s_b+1);

% Zeeman frequencies
omega_a=spin(spins{1})*B0; omega_b=spin(spins{2})*B0;

% Textbook equation, R1 first spin, rank 1 component
r1(1)=(1/6)*Ssq_b*LSqHFC*(spden(1,r_dif_c,omega_a)+...
                          spden(1,r_dif_c,omega_a-omega_b));

% Textbook equation, R1 first spin, rank 2 component
r1(1)=r1(1)+(2/27)*Ssq_b*DSqHFC*(3*spden(2,r_dif_c,omega_a)+...
                                 6*spden(2,r_dif_c,omega_a+omega_b)+...
                                 1*spden(2,r_dif_c,omega_a-omega_b));

% Textbook equation, R1 second spin, rank 1 component
r1(2)=(1/6)*Ssq_a*LSqHFC*(spden(1,r_dif_c,omega_b)+...
                          spden(1,r_dif_c,omega_a-omega_b));
 
% Textbook equation, R1 second spin, rank 2 component
r1(2)=r1(2)+(2/27)*Ssq_a*DSqHFC*(3*spden(2,r_dif_c,omega_b)+...
                                 6*spden(2,r_dif_c,omega_a+omega_b)+...
                                 1*spden(2,r_dif_c,omega_a-omega_b));

% Textbook equation, R2 first spin, rank 1 component
r2(1)=(1/12)*Ssq_b*LSqHFC*(1*spden(1,r_dif_c,omega_a)+...
                           2*spden(1,r_dif_c,omega_b)+...
                           1*spden(1,r_dif_c,omega_a-omega_b));

% Textbook equation, R2 first spin, rank 2 component
r2(1)=r2(1)+(1/27)*Ssq_b*DSqHFC*(4*spden(2,r_dif_c,0)+...
                                 3*spden(2,r_dif_c,omega_a)+...
                                 6*spden(2,r_dif_c,omega_b)+...
                                 6*spden(2,r_dif_c,omega_a+omega_b)+...
                                 1*spden(2,r_dif_c,omega_a-omega_b));

% Textbook equation, R2 second spin, rank 1 component
r2(2)=(1/12)*Ssq_a*LSqHFC*(1*spden(1,r_dif_c,omega_b)+...
                           2*spden(1,r_dif_c,omega_a)+...
                           1*spden(1,r_dif_c,omega_a-omega_b));

% Textbook equation, R2 second spin, rank 2 component
r2(2)=r2(2)+(1/27)*Ssq_a*DSqHFC*(4*spden(2,r_dif_c,0)+...
                                 3*spden(2,r_dif_c,omega_b)+...
                                 6*spden(2,r_dif_c,omega_a)+...
                                 6*spden(2,r_dif_c,omega_a+omega_b)+...
                                 1*spden(2,r_dif_c,omega_a-omega_b));

% Textbook equation, longitudinal cross-relaxation rate, rank 1 component
rx=-(1/6)*sqrt(Ssq_a*Ssq_b)*LSqHFC*spden(1,r_dif_c,omega_a-omega_b);

% Textbook equation, longitudinal cross-relaxation rate, rank 2 component
rx=rx+(2/27)*sqrt(Ssq_a*Ssq_b)*DSqHFC*(6*spden(2,r_dif_c,omega_a+omega_b)-...
                                       1*spden(2,r_dif_c,omega_a-omega_b));
    
end

% Consistency enforcement
function grumble(B0,HFC,spins,tau_c)
if (~isnumeric(B0))||(~isreal(B0))||(~isscalar(B0))
    error('B0 must be a real number.');
end
if (~isnumeric(tau_c))||(~isreal(tau_c))||...
   (~isscalar(tau_c))||(tau_c<=0)
    error('tau_c must be a positive real number.');
end
if (~isnumeric(HFC))||(~isreal(HFC))||...
   (size(HFC,1)~=3)||(size(HFC,2)~=3)
    error('HFC must be a real 3x3 matrix.');
end
if (~iscell(spins))||(numel(spins)~=2)||...
   (~ischar(spins{1}))||(~ischar(spins{2}))
    error('spins must be a cell array containing two character strings.');
end
end

% During his time at the University of Southampton, IK had formal 
% disciplinary investigations started against him for:
%
%  1. Working in his office during the Christmas break.
%
%  2. Asking a junior clerk in charge of making the useless
%     Faculty Newsletter why he existed.
%
%  3. Discussing chemicals with chemistry students inside a
%     Chemistry Department building - somebody did not like
%     the chemicals and dropped the dime.
%
%  4. Advising a student to apply for PhD positions in Cali-
%     fornia and Oxbridge instead of Southampton.
%
%  5. Refusing, "in a sarcastic and disrespectful way", to 
%     attend the ritual dressing-down session on Item 4.
%
% The cases generated mountains of paperwork; all were eventual-
% ly dropped after no wrongdoing could be established. On Item 2,
% the proceedings were adjourned because the HR could not formu-
% late the accusation. On Item 4, the case could not proceed be-
% cause the staff member accusing IK of encouraging her student
% to move to Oxford had... moved to Oxford.

