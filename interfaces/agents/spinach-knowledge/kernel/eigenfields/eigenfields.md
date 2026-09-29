# kernel/eigenfields/eigenfields.m

- Direct source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/eigenfields/eigenfields.m
- Spinach Wiki: https://spindynamics.org/wiki/index.php?title=eigenfields.m

## Signature

`tran=eigenfields(spin_system,parameters,Hz,Hc,Hmw)`

## Purpose

For `Hc+B*Hz`, find fields `B` in the requested window where an energy-level gap matches the microwave angular frequency and the transition moment through `Hmw` is significant. The formalism selects either a Hilbert-space field sweep or a Liouville-space generalised eigenproblem.

## Inputs

- `spin_system`: Spinach spin-system structure; `spin_system.bas.formalism` selects the pathway.
- `Hz`: numeric square Hermitian field-dependent Hamiltonian (Hilbert space) or commutation superoperator (Liouville space), normalised to one Tesla.
- `Hc`: numeric square Hermitian field-independent Hamiltonian or superoperator, the same size as `Hz`, containing couplings and offsets.
- `Hmw`: observable operator (Hilbert space) or observable vector (Liouville space), without the amplitude prefactor.
- `parameters.window`: two real magnetic-field endpoints in Tesla. The implementation uses their minimum and maximum; the input order does not affect the window.
- `parameters.mw_freq`: microwave frequency in Hz; the code converts it to angular frequency as `omega=2*pi*mw_freq`.
- `parameters.orientation`: three Euler angles in radians specifying system orientation.
- `parameters.tm_tol`: relative transition-moment threshold.
- `parameters.pp_tol`: peak-position tolerance in Tesla, intended to be much smaller than a typical line width.
- `parameters.fwhm`: positive transition full width at half maximum in Tesla.
- `parameters.rspt_order`: for `zeeman-hilb`, perturbation order `1`, `2`, `3`, or `4`, or `Inf` for exact diagonalisation; zero and other finite orders are rejected by `rspt_eig`. The selected order controls treatment of the off-diagonal Hamiltonian part.

The source checks that `Hc` and `Hz` are same-sized square Hermitian numeric matrices. It requires the microwave frequency, window, transition-moment tolerance, and FWHM fields; these parameter values are checked for the types and scalar/size constraints implemented in the source. It also checks peak-position tolerance as a real scalar and validates `rspt_order` for the Hilbert-space pathway. It does not impose a window endpoint ordering.

## Outputs

All transition records are filtered and sorted by increasing field.

- `tran.tf`: transition fields in Tesla, one value per retained transition.
- `tran.tm`: transition moments.
- `tran.tw`: transition widths in Tesla; the implementation assigns the supplied `parameters.fwhm`.
- `tran.pd`: energy-level population differences in the Hilbert pathway; the Liouville pathway currently assigns ones.
- `tran.ti`: transition identities. Hilbert-space rows contain source level, destination level, and that pair's branch number; Liouville-space identities are generalised-eigenvector ordinals.
- `tran.tj`: scaled field-sweep Jacobians, computed from the transition-frequency slope relative to the electron gyromagnetic factor.
- `tran.xyz`: a three-element column of `NaN` placeholders.

## Numerical / algorithmic content

For `zeeman-hilb`, the field window is initially sampled at its endpoints and its one-third and two-third points. At each knot, the routine obtains energies, eigenvectors, field derivatives, transition moments, and level populations. It tests cubic Hermite energy predictions at the one-third and two-third points against calculated, sorted energies. Intervals are trisected until both midpoint errors are below `abs(spin('E')*parameters.pp_tol)`; construction stops with an error if the knot count exceeds 1000. Eigenstates are tracked between knots by greedily assigning the largest remaining eigenvector overlaps.

A level pair is considered active if its transition moment exceeds `tm_tol` at any knot. On each interval the code forms Hermite cubic energy curves and solves their gap-minus-`omega` polynomial with [`cubic_roots.m`](cubic_roots.md). It also tests stationary points of the gap polynomial for near-zero (tangent) crossings. Transition moments, populations, and the frequency slope are interpolated at candidate roots. If either tracked state's squared overlap across the interval is not greater than 0.5, the routine rediagonalises at that root and reassigns the eigenstates by overlap before evaluating those quantities. The scaled Jacobian is `abs(spin('E'))/abs(d(energy gap)/dB)`, with `Inf` for a zero slope. The Hilbert pathway records each pair's branch number and assigns the supplied FWHM.

For `zeeman-liouv` and `sphten-liouv`, the code solves `(omega*I-Hc)*uv = B*Hz*uv`. It discards non-finite field solutions, sufficiently complex fields, generalised eigenvectors with norm below `sqrt(eps)`, solutions whose pencil residual exceeds the source's relative threshold `1e-8`, and fields outside the window. It normalises retained eigenvectors, computes transition moments as `abs(Hmw'*uv).^2`, uses the expectation of `Hz` for the field-sweep slope, and sets population differences to one.

In both pathways, transitions with moment below `parameters.tm_tol` are removed. The code's explicit consistency checks cover Hamiltonian/superoperator sizes and Hermiticity, required parameter fields, real scalar tolerances, positive FWHM, two real window values, and the Hilbert-space perturbation order. An unsupported formalism raises an error.
