# kernel/utilities/rlx_phonon.m

- Signature: `R=rlx_phonon(spin_system,H,X,I0,alpha,T,form)`

## Purpose

Constructs the spin-phonon relaxation operator in the generalised Lindblad form of Saito, Miyashita, and De Raedt (Phys. Rev. B 60, 14553 (1999)), used by Nakano and Miyashita (J. Phys. Soc. Jpn. 70, 2151 (2001)) for molecular-magnet magnetisation dynamics. The bath spectral density is `I(w)=I0*w^alpha*theta(w)`; a Hermitian coupling operator `X` gives the dissipator `d(rho)/dt=-pi*([X,R*rho]+[X,R*rho]')`. In the Hamiltonian eigenbasis, with transition frequencies `w_kn=E_k-E_n`, the thermally dressed operator has elements `R_kn=X_kn*(I(w_kn)-I(-w_kn))/(exp(hbar*w_kn/(k_B*T))-1)`.

## Physical / mathematical content

The model uses the generalised Redfield/Lindblad phonon-bath construction with a power-law spectral density. The coupling constant squared, `lambda^2`, is absorbed into `I0`, so `lambda^2*I(w)=I0*w^alpha` for positive frequency. With Hermitian `X`, the resulting Liouville-space superoperator conserves trace; its stationary relaxation destination is the thermal equilibrium state of `H` at temperature `T`, not the unit state.

## Numerical / algorithmic content

The Hamiltonian is diagonalised unless it is already diagonal; in that case the supplied basis is used directly. The code transforms and symmetrises `X` in the eigenbasis, evaluates the transition-frequency thermal spectral factor, and transforms the dressed operator back to the input basis. For `|hbar*w/(k_B*T)|<1e-3` it uses the analytic small-exponent limit; exponents above 700 are treated as infinite. In Liouville form it constructs the superoperator for the column-stretched density matrix. Because the dissipator depends on `H`, it must be rebuilt when the field changes; `pulsed_field.m` does this along the field profile and can supply diagonal `H` in its eigenbasis, avoiding a repeated diagonalisation.

## Parameters / inputs

- `spin_system` - Spinach system structure providing constants such as `hbar` and `kbol`.
- `H` - Hermitian Hilbert-space Hamiltonian in rad/s at the current magnetic field and orientation.
- `X` - Hermitian Hilbert-space spin-phonon coupling operator; `lambda` is absorbed into `I0`.
- `I0` - non-negative spectral-density prefactor such that `lambda^2*I(w)=I0*w^alpha`; units are `(rad/s)^(1-alpha)`.
- `alpha` - spectral-density exponent, 1 (Ohmic) or greater (super-Ohmic); values below 1 are unsupported because the zero-frequency limit diverges.
- `T` - positive phonon-bath temperature in kelvin.
- `form` - `'hilb'` returns the dressed coupling operator; `'liouv'` returns the relaxation superoperator.

## Outputs

- `R` - for `'hilb'`, the dressed coupling operator in the basis of `H` and `X`, used in `-pi*([X,R*rho]+[X,R*rho]')`; for `'liouv'`, the relaxation superoperator in the corresponding Liouville space, to be added as `L=H_comm+1i*R`.

## Implementation structure

Inputs are checked for Hermitian finite square matrices of matching dimension, non-negative `I0`, `alpha>=1`, positive `T`, and `form` equal to `'hilb'` or `'liouv'`. The spectral-density difference and thermal denominator are evaluated elementwise; `expm1` is used away from the small-exponent limit, and large Boltzmann exponents avoid overflow. In Liouville form the returned matrix uses the column-stretched density-matrix convention.

## References

- Saito, Miyashita, and De Raedt, *Physical Review B* **60**, 14553 (1999).
- Nakano and Miyashita, *Journal of the Physical Society of Japan* **70**, 2151 (2001).
- [Spin Dynamics Wiki: rlx_phonon.m](https://spindynamics.org/wiki/index.php?title=rlx_phonon.m)
