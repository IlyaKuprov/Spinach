# kernel/utilities/wigner_fock.m

- Signature: `W=wigner_fock(rho,alpha)`

## Purpose

Wigner function of a bosonic mode state given as a density matrix in a truncated Fock basis, evaluated at the specified points of the phase space through the displaced parity operator (Royer, Phys. Rev. A 15, 449, 1977): W(alpha)=(2/pi)*trace(rho*D(alpha)*P*D(alpha)') where D(alpha)=expm(alpha*a'-conj(alpha)*a) is the displacement operator and P=expm(1i*pi*a'*a) is the photon number parity ope- rator, both built from the ladder operators truncated to the di- mension of the density matrix. Syntax: W=wigner_fock(rho,alpha)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Each requested phase-space point costs one `expm` of the displacement generator alpha*a'-conj(alpha)*a in the truncated Fock basis and one trace against the parity operator; there is no averaging, propagation, or quadrature, and the points are evaluated independently in a loop over `numel(alpha)`.
- The only eigenvalue call is the positive-semidefiniteness test in the grumbler; the Hermiticity, unit-trace, and positivity tolerances are sqrt(eps) of the class of `rho`, so single-precision density matrices are accepted.

## Parameters / inputs

- rho -density matrix of the mode in the Fock basis
- with the levels in ascending order, [n x n]
- alpha -complex phase space coordinates, an array of
- any size; the real part is the position qua-
- drature and the imaginary part is the momen-
- tum quadrature in the units where a coherent
- state |beta> has W=(2/pi)*exp(-2*|alpha-beta|^2)

## Outputs

- W -Wigner function values, a real array of the
- same size as alpha, normalised to a unit in-
- tegral over the complex plane
- Note: the Fock basis must be large enough for the displaced
- states to fit: with the highest occupied level m of rho,
- n must be well above (sqrt(m)+|alpha|)^2 at every point,
- otherwise the truncated displacement operator distorts
- the values; pad the density matrix with zero rows and
- columns when the state itself lives in a smaller space,
- and check the unit integral on the grid in use.

## Implementation structure

- Wigner function of a bosonic mode state given as a density matrix
- in a truncated Fock basis, evaluated at the specified points of the
- phase space through the displaced parity operator (Royer, Phys.
- Rev. A 15, 449, 1977):
- W(alpha)=(2/pi)*trace(rho*D(alpha)*P*D(alpha)')
- where D(alpha)=expm(alpha*a'-conj(alpha)*a) is the displacement
- operator and P=expm(1i*pi*a'*a) is the photon number parity ope-
- rator, both built from the ladder operators truncated to the di-
- mension of the density matrix. Syntax:
- W=wigner_fock(rho,alpha)
- rho -density matrix of the mode in the Fock basis
- with the levels in ascending order, [n x n]
