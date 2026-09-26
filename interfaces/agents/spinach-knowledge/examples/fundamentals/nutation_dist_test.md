# examples/fundamentals/nutation_dist_test.m

- Signature: `nutation_dist_test()`

## Purpose

Recovery of an RF field distribution from a nutation curve measured with the same coil used for excitation and detection. A proton en- semble is driven on-resonance by a bimodal distribution of RF field amplitudes; following the reciprocity principle, the detected sig- nal of every ensemble member is weighted by its own RF field ampli- tude. The transverse magnetisation components are combined into a complex nutation

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.

## Implementation structure

- Recovery of an RF field distribution from a nutation curve measured
- with the same coil used for excitation and detection. A proton en-
- semble is driven on-resonance by a bimodal distribution of RF field
- amplitudes; following the reciprocity principle, the detected sig-
- nal of every ensemble member is weighted by its own RF field ampli-
- tude. The transverse magnetisation components are combined into a
- complex nutation curve, a receiver phase is applied, noise is add-
- ed, and nutation_dist.m is called with a user-specified Tikhonov
- regularisation parameter to recover the true nutation frequency
- distribution with the reception weight divided out.
- Calculation time: seconds
- Isotopes
