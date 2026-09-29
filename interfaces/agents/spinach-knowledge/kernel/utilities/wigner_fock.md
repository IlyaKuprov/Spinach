# kernel/utilities/wigner_fock.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/wigner_fock.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/wigner_fock.m)

## Purpose

Computes the Wigner function of a bosonic mode state given as a density matrix in a truncated Fock basis, evaluated at specified points of phase space through the displaced parity operator (Royer, Phys. Rev. A 15, 449, 1977):

```
W(alpha) = (2/pi) * trace(rho * D(alpha) * P * D(alpha)')
```

where `D(alpha) = expm(alpha*a' - conj(alpha)*a)` is the displacement operator and `P = expm(1i*pi*a'*a)` is the photon number parity operator, both built from ladder operators truncated to the dimension of the density matrix.

## Behaviour

- Syntax: `W = wigner_fock(rho, alpha)`.
- The annihilation operator is constructed as `diag(sqrt(1:(nlevels-1)), 1)` and the parity operator as `diag((-1).^(0:(nlevels-1)))`, where `nlevels` is the dimension of `rho`.
- For each point in `alpha`, the displacement operator is built via `expm(alpha(n)*an_op' - conj(alpha(n))*an_op)` and the Wigner value is `(2/pi) * real(trace(rho * disp_op * parity * disp_op'))`.
- The output `W` is a real array of the same size as `alpha`, normalised to a unit integral over the complex plane.
- Input validation (via the internal `grumble` function) requires:
  - `rho` to be a square floating-point matrix of dimension at least 2 with finite elements.
  - `rho` to be Hermitian (checked with tolerance `sqrt(eps(class(rho)))` on the 1-norm of `rho - rho'`).
  - `rho` to have unit trace (within the same tolerance).
  - `rho` to be positive semidefinite (eigenvalues of `(rho+rho')/2` not below `-tol`).
  - `alpha` to be a floating-point array with finite elements.
- The real part of the position quadrature corresponds to `real(alpha)` and the momentum quadrature to `imag(alpha)`, in units where a coherent state `|beta>` has `W = (2/pi)*exp(-2*|alpha-beta|^2)`.
- The Fock basis must be large enough for the displaced states to fit: with the highest occupied level `m` of `rho`, `n` must be well above `(sqrt(m)+|alpha|)^2` at every point, otherwise the truncated displacement operator distorts the values. The density matrix should be padded with zero rows and columns when the state itself lives in a smaller space, and the unit integral should be checked on the grid in use.

## Inputs and outputs

**Inputs**

- `rho` — density matrix of the mode in the Fock basis with the levels in ascending order, `[n x n]`.
- `alpha` — complex phase space coordinates, an array of any size; the real part is the position quadrature and the imaginary part is the momentum quadrature.

**Outputs**

- `W` — Wigner function values, a real array of the same size as `alpha`, normalised to a unit integral over the complex plane.

## References

- Royer, A. Phys. Rev. A 15, 449 (1977).
- Spin Dynamics Wiki page for `wigner_fock.m`: [https://spindynamics.org/wiki/index.php?title=wigner_fock.m](https://spindynamics.org/wiki/index.php?title=wigner_fock.m)
