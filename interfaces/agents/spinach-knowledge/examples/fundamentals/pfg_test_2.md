# examples/fundamentals/pfg_test_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/pfg_test_2.m`
- Signature: `pfg_test_2()`
- Total lines: 86

## Purpose

Demonstrate the use of the auxiliary matrix algorithm in generating a gradient sandwich multiple-quantum filter. For further details see: Calculation time: seconds

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Demonstrate the use of the auxiliary matrix algorithm in generating a
- gradient sandwich multiple-quantum filter. For further details see:
- Calculation time: seconds
- Magnet and isotopes
- Random chemical shifts and couplings
- Set the basis
- Run Spinach housekeeping
- Build the Hamiltonian
- Propagator for a pi/2 pulse
- Build initial state vector
- Determine projection quantum numbers of the basis
- Determine the coherence order of each state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `hamiltonian()`, `operator()`, `propagator()`, `rho()`, `lin2lm()`, `weighting()`, `durations_temp()`, `time_axis()`, `rho_stack()`, `grad_sandw()`, `kfigure()`, `trajan()`, `set()`.
