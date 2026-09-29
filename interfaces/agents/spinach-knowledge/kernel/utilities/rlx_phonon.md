# kernel/utilities/rlx_phonon.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rlx_phonon.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rlx_phonon.m)

## Purpose

Computes spin-phonon relaxation in the generalised Lindblad form of Saito, Miyashita, and De Raedt (Phys. Rev. B 60, 14553 (1999)), as used by Nakano and Miyashita (J. Phys. Soc. Jpn. 70, 2151 (2001)) for the magnetisation dynamics of molecular magnets. A phonon bath with spectral density `I(w)=I0*w^alpha*theta(w)` couples to the spin system through a Hermitian operator `X`.

## Behaviour

- The dissipator is `d(rho)/dt = -pi*([X,R*rho]+[X,R*rho'])`, where the thermally dressed coupling operator `R` is built in the eigenbasis of the current Hamiltonian from the transition frequencies `w_kn=(E_k-E_n)`:
  `<k|R|n> = <k|X|n> * (I(w_kn)-I(-w_kn))/(exp(hbar*w_kn/kT)-1)`.
- The square of the coupling constant `lambda` of the original papers is absorbed into the spectral density prefactor `I0`.
- If `H` is diagonal, its diagonal entries are used as eigenvalues, `X` is taken as already supplied in that eigenbasis, and diagonalisation is skipped; otherwise `H` is diagonalised with `eig` and `X` is transformed into the eigenbasis and Hermitised as `(XE+XE')/2`.
- Transition frequencies are `w=E-E.'` and Boltzmann exponents are `beta_w=hbar*w/(kbol*T)`.
- The thermal spectral density difference uses `num=I0*(max(w,0).^alpha-max(-w,0).^alpha)`; for `abs(beta_w)>=1e-3` and `beta_w<=700`, `phi=num./expm1(beta_w)`.
- Small-exponent limit (`abs(beta_w)<1e-3`) is taken analytically: `phi=I0*abs(w).^(alpha-1)*(kbol*T/hbar)-I0*sign(w).*abs(w).^alpha/2`.
- Boltzmann exponents above 700 are treated as infinite (the corresponding `phi` entries remain zero).
- The dressed operator is returned in the basis in which `H` and `X` are supplied: `R=V*(XE.*phi)*V'`.
- For `'liouv'`, the Liouville space dissipator is built for the column-stretched density matrix convention as `R=-pi*(kron(unit,X*R)-kron(X.',R)+kron((R'*X).',unit)-kron(conj(R),X))`, to be added to the Liouvillian as `L=H_comm+1i*R`.
- The dissipator depends on the Hamiltonian and must be rebuilt whenever the field changes; `pulsed_field.m` does this at every stair of the field profile, supplying `H` as the diagonal matrix of its eigenvalues and `X` in the same eigenbasis.
- The trace is conserved (the unit state is a left null vector of the superoperator) because `X` is Hermitian; the unit state itself is not stationary — the relaxation destination is the thermal equilibrium state of `H` at temperature `T`.
- Input validation (via the internal `grumble` function) requires: `H` a Hermitian matrix with finite elements; `X` a Hermitian matrix of the same dimension as `H` with finite elements; `I0` a non-negative real scalar; `alpha` a real scalar not smaller than 1 (sub-Ohmic baths make the zero-frequency limit diverge and are not supported); `T` a positive real scalar; `form` either `'hilb'` or `'liouv'`.

## Inputs and outputs

Syntax: `R=rlx_phonon(spin_system,H,X,I0,alpha,T,form)`

**Inputs**

- `spin_system` — Spinach system object providing `spin_system.tols.hbar` and `spin_system.tols.kbol`.
- `H` — Hilbert space Hamiltonian, rad/s, at the current magnetic field and orientation.
- `X` — Hilbert space spin-phonon coupling operator; the coupling constant `lambda` is absorbed into `I0`.
- `I0` — phonon spectral density prefactor times `lambda^2`, such that `lambda^2*I(w)=I0*w^alpha`; units are `(rad/s)^(1-alpha)`.
- `alpha` — spectral density exponent, 1 (Ohmic) or above (super-Ohmic).
- `T` — phonon bath temperature, Kelvin.
- `form` — `'hilb'` returns the dressed coupling operator `R`; `'liouv'` returns the relaxation superoperator.

**Outputs**

- `R` — for `'hilb'`, the dressed coupling operator in the basis in which `H` and `X` are given, such that the dissipator is `-pi*([X,R*rho]+[X,R*rho]')`; for `'liouv'`, the relaxation superoperator in the Liouville space of that basis, to be added to the Liouvillian as `L=H_comm+1i*R`.

## References

- Saito, Miyashita, De Raedt, Phys. Rev. B 60, 14553 (1999).
- Nakano and Miyashita, J. Phys. Soc. Jpn. 70, 2151 (2001).
- Spinach Wiki: [https://spindynamics.org/wiki/index.php?title=rlx_phonon.m](https://spindynamics.org/wiki/index.php?title=rlx_phonon.m)
- Source: [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rlx_phonon.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rlx_phonon.m)
