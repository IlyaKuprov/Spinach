% Solomon-Bloembergen-Morgan nuclear relaxation rates due to a
% paramagnetic centre. Syntax:
%
%       [r1,r2]=rlx_sbm(B0,nucleus,dist,a_iso,e_spin,g_eff,t1e,t2e,tau_r)
%
% Parameters:
%
%    B0       - magnet field, Tesla
%
%    nucleus  - nuclear isotope, e.g. '1H' or '13C'
%
%    dist     - electron-nucleus distance, Angstrom
%
%    a_iso    - isotropic hyperfine coupling, rad/s
%
%    e_spin   - effective electron spin quantum number
%
%    g_eff    - effective electron g-factor
%
%    t1e      - longitudinal electron relaxation time, seconds
%
%    t2e      - transverse electron relaxation time, seconds
%
%    tau_r    - rotational correlation time, seconds
%
% Outputs:
%
%    r1       - longitudinal rates [dipolar contact], Hz
%
%    r2       - transverse rates [dipolar contact], Hz
%
% The spectral density convention is J(omega,tau)=tau/(1+omega^2*tau^2).
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=rlx_sbm.m>

function [r1,r2]=rlx_sbm(B0,nucleus,dist,a_iso,e_spin,g_eff,t1e,t2e,tau_r)

% Check consistency
grumble(B0,nucleus,dist,a_iso,e_spin,g_eff,t1e,t2e,tau_r);

% Physical constants
mu0=1.25663706127e-6;
muB=9.2740100657e-24;
hbar=6.62607015e-34/(2*pi);

% Larmor frequencies
omega_i=abs(spin(nucleus)*B0);
omega_s=abs(g_eff*muB*B0/hbar);

% Effective dipolar correlation times
tau_c1=1/(1/tau_r+1/t1e);
tau_c2=1/(1/tau_r+1/t2e);

% Dipolar relaxation prefactor
spin_sq=e_spin*(e_spin+1);
dist=dist*1e-10;
pref_dd=(mu0/(4*pi))^2*(spin(nucleus)*g_eff*muB)^2*spin_sq/dist^6;

% Dipolar longitudinal relaxation rate
r1_dd=(2/15)*pref_dd*(3*tau_c1/(1+(omega_i*tau_c1)^2)+...
                      6*tau_c2/(1+((omega_i+omega_s)*tau_c2)^2)+...
                        tau_c2/(1+((omega_i-omega_s)*tau_c2)^2));

% Dipolar transverse relaxation rate
r2_dd=(1/15)*pref_dd*(4*tau_c1+...
                      3*tau_c1/(1+(omega_i*tau_c1)^2)+...
                      6*tau_c2/(1+(omega_s*tau_c2)^2)+...
                      6*tau_c2/(1+((omega_i+omega_s)*tau_c2)^2)+...
                        tau_c2/(1+((omega_i-omega_s)*tau_c2)^2));

% Contact longitudinal relaxation rate
r1_sc=(2/3)*a_iso^2*spin_sq*t2e/...
      (1+((omega_i-omega_s)*t2e)^2);

% Contact transverse relaxation rate
r2_sc=(1/3)*a_iso^2*spin_sq*(t1e+t2e/...
      (1+((omega_i-omega_s)*t2e)^2));

% Package contributions by mechanism
r1=[r1_dd r1_sc];
r2=[r2_dd r2_sc];

end

% Consistency enforcement
function grumble(B0,nucleus,dist,a_iso,e_spin,g_eff,t1e,t2e,tau_r)
if (~isnumeric(B0))||(~isreal(B0))||(~isscalar(B0))||...
   (~isfinite(B0))||(B0<=0)
    error('B0 must be a positive real finite scalar.');
end
if (~ischar(nucleus))||isempty(regexp(nucleus,'^\d+[A-Z][a-z]?\d*$','once'))
    error('nucleus must be a character string specifying a nuclear isotope.');
end
spin(nucleus);
if (~isnumeric(dist))||(~isreal(dist))||(~isscalar(dist))||...
   (~isfinite(dist))||(dist<=0)
    error('dist must be a positive real finite scalar.');
end
if (~isnumeric(a_iso))||(~isreal(a_iso))||(~isscalar(a_iso))||...
   (~isfinite(a_iso))
    error('a_iso must be a real finite scalar.');
end
if (~isnumeric(e_spin))||(~isreal(e_spin))||(~isscalar(e_spin))||...
   (~isfinite(e_spin))||(e_spin<=0)||(mod(2*e_spin,1)~=0)
    error('e_spin must be a positive integer or half-integer scalar.');
end
if (~isnumeric(g_eff))||(~isreal(g_eff))||(~isscalar(g_eff))||...
   (~isfinite(g_eff))||(g_eff<=0)
    error('g_eff must be a positive real finite scalar.');
end
if (~isnumeric(t1e))||(~isreal(t1e))||(~isscalar(t1e))||...
   (~isfinite(t1e))||(t1e<=0)
    error('t1e must be a positive real finite scalar.');
end
if (~isnumeric(t2e))||(~isreal(t2e))||(~isscalar(t2e))||...
   (~isfinite(t2e))||(t2e<=0)
    error('t2e must be a positive real finite scalar.');
end
if (~isnumeric(tau_r))||(~isreal(tau_r))||(~isscalar(tau_r))||...
   (~isfinite(tau_r))||(tau_r<=0)
    error('tau_r must be a positive real finite scalar.');
end
end

% Any sufficiently advanced incompetence is indistinguishable from malice.
%
% Grey's law


