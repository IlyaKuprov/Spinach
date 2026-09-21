% Spin-phonon relaxation superoperator in the generalised Lindblad form
% of Saito, Miyashita, and De Raedt (Phys. Rev. B 60, 14553 (1999)), as
% used by Nakano and Miyashita (J. Phys. Soc. Jpn. 70, 2151 (2001)) for
% the magnetisation dynamics of molecular magnets. A phonon bath with
% the spectral density I(w)=I0*w^alpha*theta(w) couples to the spin
% system through a Hermitian operator X; the dissipator is
%
%             d(rho)/dt = -(lambda^2*pi)*([X,R*rho]+[X,R*rho]')
%
% where R is built in the eigenbasis of the current Hamiltonian from
% the transition frequencies w_kn=(E_k-E_n):
%
%       <k|R|n> = <k|X|n> * (I(w_kn)-I(-w_kn))/(exp(hbar*w_kn/kT)-1)
%
% The Hermitian conjugate term is linear in rho for Hermitian rho, and
% so the dissipator is returned as an ordinary Liouville space super-
% operator acting on the column-stretched density matrix. Syntax:
%
%               R=rlx_phonon(spin_system,H,X,I0,alpha,T)
%
% Parameters:
%
%    H     - Hilbert space Hamiltonian, rad/s, at the current
%            magnetic field and orientation
%
%    X     - Hilbert space spin-phonon coupling operator; the
%            coupling constant lambda is absorbed into I0
%
%    I0    - phonon spectral density prefactor times lambda^2,
%            such that lambda^2*I(w)=I0*w^alpha; the units are
%            (rad/s)^(1-alpha)
%
%    alpha - spectral density exponent (sub-Ohmic below 1,
%            Ohmic at 1, super-Ohmic above 1)
%
%    T     - phonon bath temperature, Kelvin
%
% Outputs:
%
%    R     - relaxation superoperator in the Liouville space
%            of the Hilbert space in which H and X are given,
%            to be added to the Liouvillian as L=H_comm+1i*R
%
% Note: the superoperator depends on the Hamiltonian and must be
%       rebuilt whenever the field changes; pulsed_field.m does
%       this at every stair of the field profile.
%
% Note: the unit state is not damped and the trace is conserved
%       because X is Hermitian; the relaxation destination is the
%       thermal equilibrium state of H at temperature T.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=rlx_phonon.m>

function R=rlx_phonon(spin_system,H,X,I0,alpha,T)

% Check consistency
grumble(H,X,I0,alpha,T);

% Diagonalise the Hamiltonian
[V,E]=eig(full((H+H')/2),'vector');

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

% Coupling operator in the eigenbasis and the R operator back in the original basis
XE=V'*X*V; RH=V*(XE.*phi)*V';

% Liouville space dissipator, column-stretched density matrix convention
unit=speye(size(H,1));
R=-pi*(kron(unit,X*RH)-kron(X.',RH)+kron((RH'*X).',unit)-kron(conj(RH),X));

end

% Consistency enforcement
function grumble(H,X,I0,alpha,T)
if (~isnumeric(H))||(size(H,1)~=size(H,2))||any(~isfinite(H(:)))
    error('H must be a square matrix with finite elements.');
end
if (~isnumeric(X))||(~ishermitian(X))||any(size(X)~=size(H))
    error('X must be a Hermitian matrix of the same dimension as H.');
end
if (~isnumeric(I0))||(~isreal(I0))||(~isscalar(I0))||(~isfinite(I0))||(I0<0)
    error('I0 must be a non-negative real scalar.');
end
if (~isnumeric(alpha))||(~isreal(alpha))||(~isscalar(alpha))||(~isfinite(alpha))||(alpha<=0)
    error('alpha must be a positive real scalar.');
end
if (~isnumeric(T))||(~isreal(T))||(~isscalar(T))||(~isfinite(T))||(T<=0)
    error('T must be a positive real scalar.');
end
end

% I have no special talents. I am only passionately curious.
%
% Albert Einstein

