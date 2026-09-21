# kernel/utilities/phonon_oper.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/phonon_oper.m`
- Signature: `R=phonon_oper(spin_system,E,X,I0,alpha,T)`
- Total lines: 100

## Purpose

Thermally dressed spin-phonon coupling operator of the generalised Lindblad dissipator of Saito, Miyashita, and De Raedt (Phys. Rev. B 60, 14553 (1999)), in the eigenbasis of the spin Hamiltonian. For a phonon bath with the spectral density I(w)=I0*w^alpha*theta(w) that couples to the spin system through a Hermitian operator X, the dissipator is `d(rho)/dt=-pi*([X,R*rho]+[X,R*rho]')`, where R is built from the transition frequencies w_kn=(E_k-E_n) as `<k|R|n>=<k|X|n>*(I(w_kn)-I(-w_kn))/(exp(hbar*w_kn/kT)-1)`, and the square of the coupling constant lambda of the original papers is absorbed into the prefactor I0 of the spectral density. This function returns R; the corresponding Liouville space superoperator is assembled by rlx_phonon.m, and pulsed_field.m applies the dissipator as matrix products.

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- E -column vector of the eigenvalues of the Hamilto-
- nian, rad/s
- X -spin-phonon coupling operator in the eigenbasis of
- the Hamiltonian; the coupling constant lambda is
- absorbed into I0
- I0 -phonon spectral density prefactor times lambda^2,
- such that lambda^2*I(w)=I0*w^alpha; the units are
- (rad/s)^(1-alpha)
- alpha -spectral density exponent, 1 (Ohmic) or above
- (super-Ohmic); sub-Ohmic baths make the zero
- frequency limit diverge and are not supported
- T -phonon bath temperature, Kelvin

## Outputs

- R -the dressed coupling operator in the eigenbasis of
- the Hamiltonian
- Note: the thermal factor has a finite limit at zero frequency for
- alpha>=1, which is taken analytically when hbar*w/kT is below
- 1e-3; Boltzmann exponents above 700 are treated as infinite.

## Implementation structure

- Thermally dressed spin-phonon coupling operator of the generalised
- Lindblad dissipator of Saito, Miyashita, and De Raedt (Phys. Rev. B
- 60, 14553 (1999)), in the eigenbasis of the spin Hamiltonian. For a
- phonon bath with the spectral density I(w)=I0*w^alpha*theta(w) that
- couples to the spin system through a Hermitian operator X, the dis-
- sipator is
- d(rho)/dt = -pi*([X,R*rho]+[X,R*rho]')
- where R is built from the transition frequencies w_kn=(E_k-E_n):
- <k|R|n> = <k|X|n> * (I(w_kn)-I(-w_kn))/(exp(hbar*w_kn/kT)-1)
- and the square of the coupling constant lambda of the original
- papers is absorbed into the prefactor I0 of the spectral density.
- This function returns R; the corresponding Liouville space super-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `phi()`, `num()`, `expm1()`, `beta_w()`, `sign()`, `iscolumn()`, `any()`, `ishermitian()`, `isscalar()`.
