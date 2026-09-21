% Spin-phonon relaxation in the generalised Lindblad form of Saito,
% Miyashita, and De Raedt (Phys. Rev. B 60, 14553 (1999)), as used by
% Nakano and Miyashita (J. Phys. Soc. Jpn. 70, 2151 (2001)) for the
% magnetisation dynamics of molecular magnets. A phonon bath with the
% spectral density I(w)=I0*w^alpha*theta(w) couples to the spin system
% through a Hermitian operator X; the dissipator is
%
%                d(rho)/dt = -pi*([X,R*rho]+[X,R*rho]')
%
% where the thermally dressed coupling operator R is built in the ei-
% genbasis of the current Hamiltonian from the transition frequencies
% w_kn=(E_k-E_n):
%
%       <k|R|n> = <k|X|n> * (I(w_kn)-I(-w_kn))/(exp(hbar*w_kn/kT)-1)
%
% and the square of the coupling constant lambda of the original
% papers is absorbed into the prefactor I0 of the spectral density.
%
% The function returns either the dressed operator R in the basis in
% which H and X are supplied, for dissipators applied as Hilbert space
% matrix products, or the Liouville space superoperator of the dissi-
% pator: the Hermitian conjugate term is linear in rho for Hermitian
% rho, and so the dissipator is an ordinary superoperator acting on
% the column-stretched density matrix. Syntax:
%
%             R=rlx_phonon(spin_system,H,X,I0,alpha,T,form)
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
%    alpha - spectral density exponent, 1 (Ohmic) or above
%            (super-Ohmic); sub-Ohmic baths make the zero
%            frequency limit diverge and are not supported
%
%    T     - phonon bath temperature, Kelvin
%
%    form  - 'hilb' returns the dressed coupling operator R,
%            'liouv' returns the relaxation superoperator
%
% Outputs:
%
%    R     - for 'hilb', the dressed coupling operator in the
%            basis in which H and X are given, such that the
%            dissipator is -pi*([X,R*rho]+[X,R*rho]'); for
%            'liouv', the relaxation superoperator in the
%            Liouville space of that basis, to be added to
%            the Liouvillian as L=H_comm+1i*R
%
% Note: the dissipator depends on the Hamiltonian and must be re-
%       built whenever the field changes; pulsed_field.m does this
%       at every stair of the field profile, supplying H as the
%       diagonal matrix of its eigenvalues and X in the same eigen-
%       basis, in which case the diagonalisation is skipped and
%       the dressed operator is returned in that eigenbasis.
%
% Note: the thermal factor has a finite limit at zero frequency for
%       alpha>=1, which is taken analytically when hbar*w/kT is below
%       1e-3; Boltzmann exponents above 700 are treated as infinite.
%
% Note: the trace is conserved (the unit state is a left null vector
%       of the superoperator) because X is Hermitian; the unit state
%       itself is not stationary, the relaxation destination is the
%       thermal equilibrium state of H at temperature T.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=rlx_phonon.m>

function R=rlx_phonon(spin_system,H,X,I0,alpha,T,form)

% Check consistency
grumble(H,X,I0,alpha,T,form);

% Eigensystem of the Hamiltonian and the coupling operator in the eigenbasis
if isdiag(H)
    V=speye(size(H)); E=full(diag(H)); XE=X;
else
    [V,E]=eig(full(H),'vector'); XE=V'*X*V; XE=(XE+XE')/2;
end

% Transition frequencies and Boltzmann exponents
w=E-E.'; beta_w=spin_system.tols.hbar*w/(spin_system.tols.kbol*T);

% Thermal spectral density difference, with the small and large exponent limits
num=I0*(max(w,0).^alpha-max(-w,0).^alpha); phi=zeros(size(w));
normal=(abs(beta_w)>=1e-3)&(beta_w<=700); small=abs(beta_w)<1e-3;
phi(normal)=num(normal)./expm1(beta_w(normal));
phi(small)=I0*abs(w(small)).^(alpha-1)*(spin_system.tols.kbol*T/spin_system.tols.hbar)-...
           I0*sign(w(small)).*abs(w(small)).^alpha/2;

% Dressed coupling operator back in the basis of the Hamiltonian
R=V*(XE.*phi)*V';

% Liouville space dissipator, column-stretched density matrix convention
if strcmp(form,'liouv')
    unit=speye(size(H,1));
    R=-pi*(kron(unit,X*R)-kron(X.',R)+kron((R'*X).',unit)-kron(conj(R),X));
end

end

% Consistency enforcement
function grumble(H,X,I0,alpha,T,form)
if (~isnumeric(H))||(~ishermitian(H))||any(~isfinite(H(:)))
    error('H must be a Hermitian matrix with finite elements.');
end
if (~isnumeric(X))||(~ishermitian(X))||any(size(X)~=size(H))||any(~isfinite(X(:)))
    error('X must be a Hermitian matrix of the same dimension as H with finite elements.');
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
if (~ischar(form))||(~ismember(form,{'hilb','liouv'}))
    error('form must be ''hilb'' or ''liouv''.');
end
end

% I have no special talents. I am only passionately curious.
%
% Albert Einstein

