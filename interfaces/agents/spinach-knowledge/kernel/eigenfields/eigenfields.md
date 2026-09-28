# kernel/eigenfields/eigenfields.m

- Signature: `tran=eigenfields(spin_system,parameters,Hz,Hc,Hmw)`

## Purpose

Computes resonance fields for `Hc+B*Hz`: fields `B` where the difference between two eigenvalues equals the supplied microwave frequency and the transition moment across `Hmw` is significant.

## Parameters / inputs

- `spin_system` — spin-system specification; its formalism selects the Hilbert- or Liouville-space pathway.
- `Hz` — field-dependent part of the laboratory-frame Hamiltonian (Hilbert space) or commutation superoperator (Liouville space), normalised to 1 Tesla.
- `Hc` — field-independent part of that operator, containing couplings and offsets.
- `Hmw` — observable operator (Hilbert space) or observable vector (Liouville space), without the amplitude prefactor.
- `parameters.window` — magnetic-field window in Tesla, given by two real elements.
- `parameters.mw_freq` — microwave frequency in Hz.
- `parameters.orientation` — three Euler angles in radians specifying the system orientation.
- `parameters.tm_tol` — relative transition-moment tolerance.
- `parameters.pp_tol` — peak-position tolerance in Tesla; it should be much smaller than the typical line width.
- `parameters.fwhm` — transition full width at half maximum in Tesla.
- `parameters.rspt_order` — perturbation-theory order for the off-diagonal part of the Hamiltonian; `Inf` requests exact diagonalisation. Required for the `zeeman-hilb` pathway.

## Outputs

- `tran.tf` — transition fields in Tesla.
- `tran.tm` — transition moments.
- `tran.tw` — transition FWHMs in Tesla.
- `tran.pd` — energy-level population differences in the Hilbert-space pathway; set to one in the Liouville-space pathway.
- `tran.ti` — transition identities: source level, destination level, and branch number in Hilbert space; generalised-eigenvector ordinal in Liouville space.
- `tran.tj` — scaled field-sweep Jacobians.
- `tran.xyz` — placeholder coordinates, set to three `NaN` values.

## Numerical / algorithmic content

- The `zeeman-hilb` pathway starts with a four-point grid spanning `parameters.window`. It repeatedly trisects unconverged intervals, comparing Hermite-spline energy predictions with eigensystem calculations until the energy errors meet a tolerance derived from `parameters.pp_tol`; grid construction fails if the grid exceeds 1,000 points. It tracks eigenstates between knots by maximum eigenvector overlap and marks level pairs active when their transition moment exceeds `parameters.tm_tol` at a grid point.
- For active pairs, the Hilbert-space pathway finds roots of cubic interpolants of the transition-frequency gap, including near-zero gaps at spline extrema. It interpolates transition moments and population differences, but recalculates them and the frequency slope at a root when state tracking across its interval is unstable. Each retained root gets a field, a scaled Jacobian, the supplied `parameters.fwhm`, and a pair-specific branch identity.
- The `zeeman-liouv` and `sphten-liouv` pathways solve the generalised eigenproblem `(omega*I-Hc)*uv = B*Hz*uv`, where `omega=2*pi*parameters.mw_freq`. They discard non-finite, appreciably complex, high-residual, and out-of-window field solutions. After normalising the eigenvectors, they calculate moments as `abs(Hmw'*uv).^2`, assign unit population differences, use `parameters.fwhm` for widths, and identify transitions by eigenvector ordinal.
- Both pathways discard moments below `parameters.tm_tol` and sort the remaining transitions by field. Inputs are checked for compatible, Hermitian `Hc` and `Hz`, required parameters, and valid numeric values; an unsupported formalism raises an error.

## Reference

- https://spindynamics.org/wiki/index.php?title=eigenfields.m