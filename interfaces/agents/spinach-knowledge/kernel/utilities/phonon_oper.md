# kernel/utilities/phonon_oper.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/phonon_oper.m`
- Signature: `R=phonon_oper(spin_system,E,X,I0,alpha,T)`
- Total lines: 96

## Purpose

Thermally dressed spin-phonon coupling operator of the generalised Lindblad dissipator of Saito, Miyashita, and De Raedt (Phys. Rev. B 60, 14553 (1999)), in the eigenbasis of the spin Hamiltonian. For a phonon bath with the spectral density I(w)=I0*w^alpha*theta(w) that couples to the spin system through a Hermitian operator X, the dis- sipator is d(rho)/dt = -(lambda^2*pi)*([X,R*rho]+[X,R*rho]') where R is built fro

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 53-54: Check consistency; implemented by `grumble(E,X,I0,alpha,T)`.
- Lines 56-57: Transition frequencies in rad/s; implemented by `w=E-E.'`.
- Lines 59-60: Boltzmann exponents; implemented by `beta_w=spin_system.tols.hbar*w/(spin_system.tols.kbol*T)`.
- Lines 62-63: Thermal spectral density difference, with the small and large exponent limits; implemented by `num=I0*(max(w,0).^alpha-max(-w,0).^alpha); phi=zeros(size(w))`.
- Lines 69-70: Dressed coupling operator; implemented by `R=X.*phi`.

### Key state/data transformations

- Lines 57: computes `w` using `w=E-E.'`.
- Lines 60: computes `beta_w` using `beta_w=spin_system.tols.hbar*w/(spin_system.tols.kbol*T)`.
- Lines 63: computes `num` using `num=I0*(max(w,0).^alpha-max(-w,0).^alpha); phi=zeros(size(w))`.
- Lines 65: computes `phi(normal)` using `phi(normal)=num(normal)./expm1(beta_w(normal))`.
- Lines 66-67: computes `phi(small)` using `phi(small)=I0*abs(w(small)).^(alpha-1)*(spin_system.tols.kbol*T/spin_system.tols.hbar)- I0*sign(w(small)).*abs(w(small)).^alpha/2`.
- Lines 70: computes `R` using `R=X.*phi`.

### Local helper functions

- Line 75: `grumble()` — `function grumble(E,X,I0,alpha,T)`.
  - Representative operation: `if (~isnumeric(E))||(~isreal(E))||(~iscolumn(E))||any(~isfinite(E))`.
  - Representative operation: `error('E must be a column vector of real finite eigenvalues.')`.

## Parameters / inputs

- E -column vector of the eigenvalues of the Hamilto-
- nian, rad/s
- X -spin-phonon coupling operator in the eigenbasis of
- the Hamiltonian; the coupling constant lambda is
- absorbed into I0
- I0 -phonon spectral density prefactor times lambda^2,
- such that lambda^2*I(w)=I0*w^alpha; the units are
- (rad/s)^(1-alpha)
- alpha -spectral density exponent (sub-Ohmic below 1,
- Ohmic at 1, super-Ohmic above 1)
- T -phonon bath temperature, Kelvin

## Outputs

- R -the dressed coupling operator in the eigenbasis of
- the Hamiltonian
- Note: the thermal factor has a finite limit at zero frequency for
- alpha>=1, which is taken analytically when hbar*w/kT is below
- 1e-3; exponents above 700 are treated as infinite.

## Implementation structure

- Thermally dressed spin-phonon coupling operator of the generalised
- Lindblad dissipator of Saito, Miyashita, and De Raedt (Phys. Rev. B
- 60, 14553 (1999)), in the eigenbasis of the spin Hamiltonian. For a
- phonon bath with the spectral density I(w)=I0*w^alpha*theta(w) that
- couples to the spin system through a Hermitian operator X, the dis-
- sipator is
- d(rho)/dt = -(lambda^2*pi)*([X,R*rho]+[X,R*rho]')
- where R is built from the transition frequencies w_kn=(E_k-E_n):
- <k|R|n> = <k|X|n> * (I(w_kn)-I(-w_kn))/(exp(hbar*w_kn/kT)-1)
- This function returns R; the corresponding Liouville space super-
- operator is assembled by rlx_phonon.m, and pulsed_field.m applies
- the dissipator as matrix products. Syntax:

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `phi()`, `num()`, `expm1()`, `beta_w()`, `sign()`, `iscolumn()`, `any()`, `ishermitian()`, `isscalar()`.
