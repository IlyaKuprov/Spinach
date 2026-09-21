# kernel/utilities/rlx_phonon.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/rlx_phonon.m`
- Signature: `R=rlx_phonon(spin_system,H,X,I0,alpha,T)`
- Total lines: 108

## Purpose

Spin-phonon relaxation superoperator in the generalised Lindblad form of Saito, Miyashita, and De Raedt (Phys. Rev. B 60, 14553 (1999)), as used by Nakano and Miyashita (J. Phys. Soc. Jpn. 70, 2151 (2001)) for the magnetisation dynamics of molecular magnets. A phonon bath with the spectral density I(w)=I0*w^alpha*theta(w) couples to the spin system through a Hermitian operator X; the dissipator is d(rho)/dt = -(lambd

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 58-59: Check consistency; implemented by `grumble(H,X,I0,alpha,T)`.
- Lines 61-62: Diagonalise the Hamiltonian; implemented by `[V,E]=eig(full((H+H')/2),'vector')`.
- Lines 64-65: Transition frequencies in rad/s; implemented by `w=E-E.'`.
- Lines 67-68: Boltzmann exponents; implemented by `beta_w=spin_system.tols.hbar*w/(spin_system.tols.kbol*T)`.
- Lines 70-71: Thermal spectral density difference, with the small and large exponent limits; implemented by `num=I0*(max(w,0).^alpha-max(-w,0).^alpha); phi=zeros(size(w))`.
- Lines 77-78: Coupling operator in the eigenbasis and the R operator back in the original basis; implemented by `XE=V'*X*V; RH=V*(XE.*phi)*V'`.
- Lines 80-81: Liouville space dissipator, column-stretched density matrix convention; implemented by `unit=speye(size(H,1))`.

### Key state/data transformations

- Lines 62: computes `[V,E]` using `[V,E]=eig(full((H+H')/2),'vector')`.
- Lines 65: computes `w` using `w=E-E.'`.
- Lines 68: computes `beta_w` using `beta_w=spin_system.tols.hbar*w/(spin_system.tols.kbol*T)`.
- Lines 71: computes `num` using `num=I0*(max(w,0).^alpha-max(-w,0).^alpha); phi=zeros(size(w))`.
- Lines 73: computes `phi(normal)` using `phi(normal)=num(normal)./expm1(beta_w(normal))`.
- Lines 74-75: computes `phi(small)` using `phi(small)=I0*abs(w(small)).^(alpha-1)*(spin_system.tols.kbol*T/spin_system.tols.hbar)- I0*sign(w(small)).*abs(w(small)).^alpha/2`.
- Lines 78: computes `XE` using `XE=V'*X*V; RH=V*(XE.*phi)*V'`.
- Lines 81: computes `unit` using `unit=speye(size(H,1))`.
- Lines 82: computes `R` using `R=-pi*(kron(unit,X*RH)-kron(X.',RH)+kron((RH'*X).',unit)-kron(conj(RH),X))`.

### Local helper functions

- Line 87: `grumble()` — `function grumble(H,X,I0,alpha,T)`.
  - Representative operation: `if (~isnumeric(H))||(size(H,1)~=size(H,2))||any(~isfinite(H(:)))`.
  - Representative operation: `error('H must be a square matrix with finite elements.')`.

## Parameters / inputs

- H -Hilbert space Hamiltonian, rad/s, at the current
- magnetic field and orientation
- X -Hilbert space spin-phonon coupling operator; the
- coupling constant lambda is absorbed into I0
- I0 -phonon spectral density prefactor times lambda^2,
- such that lambda^2*I(w)=I0*w^alpha; the units are
- (rad/s)^(1-alpha)
- alpha -spectral density exponent (sub-Ohmic below 1,
- Ohmic at 1, super-Ohmic above 1)
- T -phonon bath temperature, Kelvin

## Outputs

- R -relaxation superoperator in the Liouville space
- of the Hilbert space in which H and X are given,
- to be added to the Liouvillian as L=H_comm+1i*R
- Note: the superoperator depends on the Hamiltonian and must be
- rebuilt whenever the field changes; pulsed_field.m does
- this at every stair of the field profile.
- Note: the unit state is not damped and the trace is conserved
- because X is Hermitian; the relaxation destination is the
- thermal equilibrium state of H at temperature T.

## Implementation structure

- Spin-phonon relaxation superoperator in the generalised Lindblad form
- of Saito, Miyashita, and De Raedt (Phys. Rev. B 60, 14553 (1999)), as
- used by Nakano and Miyashita (J. Phys. Soc. Jpn. 70, 2151 (2001)) for
- the magnetisation dynamics of molecular magnets. A phonon bath with
- the spectral density I(w)=I0*w^alpha*theta(w) couples to the spin
- system through a Hermitian operator X; the dissipator is
- d(rho)/dt = -(lambda^2*pi)*([X,R*rho]+[X,R*rho]')
- where R is built in the eigenbasis of the current Hamiltonian from
- the transition frequencies w_kn=(E_k-E_n):
- <k|R|n> = <k|X|n> * (I(w_kn)-I(-w_kn))/(exp(hbar*w_kn/kT)-1)
- The Hermitian conjugate term is linear in rho for Hermitian rho, and
- so the dissipator is returned as an ordinary Liouville space super-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `phi()`, `num()`, `expm1()`, `beta_w()`, `sign()`, `speye()`, `conj()`, `any()`, `ishermitian()`, `isscalar()`.
