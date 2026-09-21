% Thermally dressed spin-phonon coupling operator of the generalised
% Lindblad dissipator of Saito, Miyashita, and De Raedt (Phys. Rev. B
% 60, 14553 (1999)), in the eigenbasis of the spin Hamiltonian. For a
% phonon bath with the spectral density I(w)=I0*w^alpha*theta(w) that
% couples to the spin system through a Hermitian operator X, the dis-
% sipator is
%
%                d(rho)/dt = -pi*([X,R*rho]+[X,R*rho]')
%
% where R is built from the transition frequencies w_kn=(E_k-E_n):
%
%       <k|R|n> = <k|X|n> * (I(w_kn)-I(-w_kn))/(exp(hbar*w_kn/kT)-1)
%
% and the square of the coupling constant lambda of the original
% papers is absorbed into the prefactor I0 of the spectral density.
%
% This function returns R; the corresponding Liouville space super-
% operator is assembled by rlx_phonon.m, and pulsed_field.m applies
% the dissipator as matrix products. Syntax:
%
%               R=phonon_oper(spin_system,E,X,I0,alpha,T)
%
% Parameters:
%
%    E     - column vector of the eigenvalues of the Hamilto-
%            nian, rad/s
%
%    X     - spin-phonon coupling operator in the eigenbasis of
%            the Hamiltonian; the coupling constant lambda is
%            absorbed into I0
%
%    I0    - phonon spectral density prefactor times lambda^2,
%            such that lambda^2*I(w)=I0*w^alpha; the units are
%            (rad/s)^(1-alpha)
%
%    alpha - spectral density exponent, 1 (Ohmic) or above
%            (super-Ohmic); sub-Ohmic baths make the zero
%            frequency limit diverge and are not supported
%
%    T     - phonon bath temperature, Kelvin
%
% Outputs:
%
%    R     - the dressed coupling operator in the eigenbasis of
%            the Hamiltonian
%
% Note: the thermal factor has a finite limit at zero frequency for
%       alpha>=1, which is taken analytically when hbar*w/kT is below
%       1e-3; Boltzmann exponents above 700 are treated as infinite.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=phonon_oper.m>

function R=phonon_oper(spin_system,E,X,I0,alpha,T)

% Check consistency
grumble(E,X,I0,alpha,T);

% Transition frequencies in rad/s
w=E-E.';

% Boltzmann exponents
beta_w=spin_system.tols.hbar*w/(spin_system.tols.kbol*T);

% Thermal spectral density difference, with the small and large exponent limits
num=I0*(max(w,0).^alpha-max(-w,0).^alpha); phi=zeros(size(w));
normal=(abs(beta_w)>=1e-3)&(beta_w<=700); small=abs(beta_w)<1e-3;
phi(normal)=num(normal)./expm1(beta_w(normal));
phi(small)=I0*abs(w(small)).^(alpha-1)*(spin_system.tols.kbol*T/spin_system.tols.hbar)-...
           I0*sign(w(small)).*abs(w(small)).^alpha/2;

% Dressed coupling operator
R=X.*phi;

end

% Consistency enforcement
function grumble(E,X,I0,alpha,T)
if (~isnumeric(E))||(~isreal(E))||(~iscolumn(E))||any(~isfinite(E))
    error('E must be a column vector of real finite eigenvalues.');
end
if (~isnumeric(X))||(~ishermitian(X))||any(size(X)~=numel(E))
    error('X must be a Hermitian matrix of the same dimension as E.');
end
if (~isnumeric(I0))||(~isreal(I0))||(~isscalar(I0))||(~isfinite(I0))||(I0<0)
    error('I0 must be a non-negative real scalar.');
end
if (~isnumeric(alpha))||(~isreal(alpha))||(~isscalar(alpha))||(~isfinite(alpha))||(alpha<1)
    error('alpha must be a real scalar not smaller than 1.');
end
if (~isnumeric(T))||(~isreal(T))||(~isscalar(T))||(~isfinite(T))||(T<=0)
    error('T must be a positive real scalar.');
end
end

% The purpose of computing is insight, not numbers.
%
% Richard Hamming

