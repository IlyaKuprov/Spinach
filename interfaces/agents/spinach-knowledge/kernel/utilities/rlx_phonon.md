# kernel/utilities/rlx_phonon.m

- Signature: `R=rlx_phonon(spin_system,H,X,I0,alpha,T,form)`

## Purpose

Spin-phonon relaxation in the generalised Lindblad form of Saito, Miyashita, and De Raedt (Phys. Rev. B 60, 14553 (1999)), as used by Nakano and Miyashita (J. Phys. Soc. Jpn. 70, 2151 (2001)) for the magnetisation dynamics of molecular magnets. A phonon bath with the spectral density I(w)=I0*w^alpha*theta(w) couples to the spin system through a Hermitian operator X; the dissipator is `d(rho)/dt=-pi*([X,R*rho]+[X,R*rho]')`, where the thermally dressed coupling operator R is built in the eigenbasis of the current Hamiltonian from the transition frequencies w_kn=(E_k-E_n) as `<k|R|n>=<k|X|n>*(I(w_kn)-I(-w_kn))/(exp(hbar*w_kn/kT)-1)`, and the square of the coupling constant lambda of the original papers is absorbed into the prefactor I0 of the spectral density. The function returns either the dressed operator R in the basis in which H and X are supplied, for dissipators applied as Hilbert space matrix products, or the Liouville space superoperator of the dissipator: the Hermitian conjugate term is linear in rho for Hermitian rho, and so the dissipator is an ordinary superoperator acting on the column-stretched density matrix.

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator. When H is already diagonal, the diagonalisation is skipped and the supplied basis is taken as the eigenbasis.
- The thermal factor `(I(w)-I(-w))/(exp(hbar*w/kT)-1)` is evaluated with `expm1` at ordinary exponents, replaced by its analytic zero-frequency expansion where `|hbar*w/kT|<1e-3`, and set to zero where the exponent exceeds 700.
- `phi`, `num`, and `beta_w` are arrays indexed by logical masks, not function calls.

## Parameters / inputs

- H -Hilbert space Hamiltonian, rad/s, at the current
- magnetic field and orientation
- X -Hilbert space spin-phonon coupling operator; the
- coupling constant lambda is absorbed into I0
- I0 -phonon spectral density prefactor times lambda^2,
- such that lambda^2*I(w)=I0*w^alpha; the units are
- (rad/s)^(1-alpha)
- alpha -spectral density exponent, 1 (Ohmic) or above
- (super-Ohmic); sub-Ohmic baths make the zero
- frequency limit diverge and are not supported
- T -phonon bath temperature, Kelvin
- form -'hilb' returns the dressed coupling operator R,
- 'liouv' returns the relaxation superoperator

## Outputs

- R -for 'hilb', the dressed coupling operator in the
- basis in which H and X are given, such that the
- dissipator is -pi*([X,R*rho]+[X,R*rho]'); for
- 'liouv', the relaxation superoperator in the
- Liouville space of that basis, to be added to
- the Liouvillian as L=H_comm+1i*R
- Note: the dissipator depends on the Hamiltonian and must be re-
- built whenever the field changes; pulsed_field.m does this
- at every stair of the field profile, supplying H as the
- diagonal matrix of its eigenvalues and X in the same eigen-
- basis, in which case the diagonalisation is skipped and
- the dressed operator is returned in that eigenbasis.
- Note: the thermal factor has a finite limit at zero frequency for
- alpha>=1, which is taken analytically when hbar*w/kT is below
- 1e-3; Boltzmann exponents above 700 are treated as infinite.
- Note: the trace is conserved (the unit state is a left null vector
- of the superoperator) because X is Hermitian; the unit state
- itself is not stationary, the relaxation destination is the
- thermal equilibrium state of H at temperature T.

## Implementation structure

- Spin-phonon relaxation in the generalised Lindblad form of Saito,
- Miyashita, and De Raedt (Phys. Rev. B 60, 14553 (1999)), as used by
- Nakano and Miyashita (J. Phys. Soc. Jpn. 70, 2151 (2001)) for the
- magnetisation dynamics of molecular magnets. A phonon bath with the
- spectral density I(w)=I0*w^alpha*theta(w) couples to the spin system
- through a Hermitian operator X; the dissipator is
- d(rho)/dt = -pi*([X,R*rho]+[X,R*rho]')
- where the thermally dressed coupling operator R is built in the ei-
- genbasis of the current Hamiltonian from the transition frequencies
- w_kn=(E_k-E_n):
- <k|R|n> = <k|X|n> * (I(w_kn)-I(-w_kn))/(exp(hbar*w_kn/kT)-1)
- and the square of the coupling constant lambda of the original
