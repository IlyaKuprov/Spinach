% Redfield theory expressions for quadrupolar relaxation rates,
% isotropic tumbling in liquid phase. Syntax:
%
%         [r1,r2,t1,t2]=rlx_nqi(I,omega,C_q,eta_q,tau_c)
%
% Parameters:
%
%    I      - nuclear spin quantum number
%
%    omega  - nuclear Zeeman frequency, rad/s
%
%    C_q    - quadrupolar coupling constant, 
%             e^2*q*Q/h in Hz
%
%    eta_q  - quadrupolar tensor asymmetry
%
%    tau_c  - rotational correlation time, seconds
%
% Outputs:
%
%    r1     - longitudinal relaxation rate, Hz
%
%    r2     - transverse relaxation rate, Hz
%
%    t1     - longitudinal relaxation time, seconds
%
%    t2     - transverse relaxation time, seconds
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rlx_nqi.m>

function [r1,r2,t1,t2]=rlx_nqi(I,omega,C_q,eta_q,tau_c)

% Check consistency
grumble(I,omega,C_q,eta_q,tau_c);

% Get the quadrupolar tensor in rad/s
Q=2*pi*eeqq2nqi(C_q,eta_q,I,[0 0 0]);

% Blicharsky invariant and diffusion coefficient
[~,DSqQ]=blinv(Q); D=1/(6*tau_c);

% Longitudinal relaxation rate
r1=(2/15)*(2*I-1)*(2*I+3)*DSqQ*(1*spden(2,D,omega)+...
                                4*spden(2,D,2*omega));

% Transverse relaxation rate
r2=(1/15)*(2*I-1)*(2*I+3)*DSqQ*(3*spden(2,D,0)+...
                                5*spden(2,D,omega)+...
                                2*spden(2,D,2*omega));

% Convert rates to times
t1=1/r1; t2=1/r2;

end

% Consistency enforcement
function grumble(I,omega,C_q,eta_q,tau_c)
if (~isnumeric(C_q))||(~isnumeric(eta_q))||(~isnumeric(I))||...
   (~isnumeric(omega))||(~isnumeric(tau_c))
    error('all inputs must be numeric.');
end
if (~isreal(C_q))||(~isreal(eta_q))||(~isreal(I))||...
   (~isreal(omega))||(~isreal(tau_c))
    error('all inputs must be real.');
end
if (~isscalar(C_q))||(~isscalar(eta_q))||(~isscalar(I))||...
   (~isscalar(omega))||(~isscalar(tau_c))
    error('all inputs must be scalars.');
end
if (numel(I)~=1)||(I<1)||(mod(2*I+1,1)~=0)
    error('I must be an integer or half-integer greater or equal to 1.');
end
if tau_c<=0
    error('tau_c must be positive.');
end
end

% The problem with today's world is that everyone believes they
% have the right to express their opinion AND have others listen
% to it. The correct statement of individual rights is that every-
% one has the right to an opinion, but crucially, that opinion 
% can be roundly ignored and even made fun of, particularly if it
% is demonstrably nonsense!
%
% Brian Cox

