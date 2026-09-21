# kernel/utilities/rlx_phonon.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/rlx_phonon.m`
- Signature: `R=rlx_phonon(spin_system,H,X,I0,alpha,T)`
- Total lines: 101

## Purpose

Spin-phonon relaxation superoperator in the generalised Lindblad form of Saito, Miyashita, and De Raedt (Phys. Rev. B 60, 14553 (1999)), as used by Nakano and Miyashita (J. Phys. Soc. Jpn. 70, 2151 (2001)) for the magnetisation dynamics of molecular magnets. A phonon bath with the spectral density I(w)=I0*w^alpha*theta(w) couples to the spin system through a Hermitian operator X; the dissipator is `d(rho)/dt=-pi*([X,R*rho]+[X,R*rho]')`, where R is built in the eigenbasis of the current Hamiltonian from the transition frequencies w_kn=(E_k-E_n) as `<k|R|n>=<k|X|n>*(I(w_kn)-I(-w_kn))/(exp(hbar*w_kn/kT)-1)`, and the square of the coupling constant lambda of the original papers is absorbed into the prefactor I0 of the spectral density. The Hermitian conjugate term is linear in rho for Hermitian rho, and so the dissipator is returned as an ordinary Liouville space superoperator acting on the column-stretched density matrix. The R operator itself is built by phonon_oper.m.

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Outputs

- R -relaxation superoperator in the Liouville space
- of the Hilbert space in which H and X are given,
- to be added to the Liouvillian as L=H_comm+1i*R
- Note: the superoperator depends on the Hamiltonian and must be
- rebuilt whenever the field changes; pulsed_field.m does
- this at every stair of the field profile.
- Note: the trace is conserved (the unit state is a left null vector
- of R) because X is Hermitian; the unit state itself is not
- stationary, the relaxation destination is the thermal equi-
- librium state of H at temperature T.

## Implementation structure

- Spin-phonon relaxation superoperator in the generalised Lindblad form
- of Saito, Miyashita, and De Raedt (Phys. Rev. B 60, 14553 (1999)), as
- used by Nakano and Miyashita (J. Phys. Soc. Jpn. 70, 2151 (2001)) for
- the magnetisation dynamics of molecular magnets. A phonon bath with
- the spectral density I(w)=I0*w^alpha*theta(w) couples to the spin
- system through a Hermitian operator X; the dissipator is
- d(rho)/dt = -pi*([X,R*rho]+[X,R*rho]')
- where R is built in the eigenbasis of the current Hamiltonian from
- the transition frequencies w_kn=(E_k-E_n):
- <k|R|n> = <k|X|n> * (I(w_kn)-I(-w_kn))/(exp(hbar*w_kn/kT)-1)
- and the square of the coupling constant lambda of the original
- papers is absorbed into the prefactor I0 of the spectral density.

## Internal Spinach / MATLAB structure cues

- Called routines in the main body: `grumble()`, `eig()`, `full()`, `phonon_oper()`, `speye()`, `kron()`, `conj()`.
