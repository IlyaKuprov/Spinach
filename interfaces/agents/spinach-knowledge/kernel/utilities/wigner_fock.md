# kernel/utilities/wigner_fock.m

- Signature: `W=wigner_fock(rho,alpha)`

## Purpose

Computes the Wigner function of a bosonic mode density matrix in a truncated Fock basis at specified phase-space points using displaced parity (Royer, Phys. Rev. A 15, 449, 1977).

## Parameters

- `rho`: An `n x n` density matrix with Fock levels in ascending order. It must be a finite floating-point matrix with `n >= 2`, and must be Hermitian, positive semidefinite, and unit trace within numerical tolerance.
- `alpha`: A finite floating-point array of complex phase-space coordinates of any size. Its real and imaginary parts are the position and momentum quadratures, respectively, in units where a coherent state `|beta>` has `W(alpha)=(2/pi)*exp(-2*|alpha-beta|^2)`.

## Output

- `W`: Real Wigner-function values with the same size as `alpha`, normalized to unit integral over the complex plane.

## Method and truncation

The function evaluates `W(alpha)=(2/pi)*trace(rho*D(alpha)*P*D(alpha)')`, where `D(alpha)=expm(alpha*a'-conj(alpha)*a)` and `P=expm(1i*pi*a'*a)` is photon-number parity. The ladder operators are truncated to the dimension of `rho`; each point is evaluated independently using a matrix exponential.

The Fock basis must accommodate the displaced state. If `m` is the highest occupied level, `n` should be well above `(sqrt(m)+|alpha|)^2` at every point. Otherwise, the truncated displacement distorts the result. Pad `rho` with zero rows and columns when needed, and check the unit integral on the grid used.

Source: <https://spindynamics.org/wiki/index.php?title=wigner_fock.m> · ilya.kuprov@weizmann.ac.il