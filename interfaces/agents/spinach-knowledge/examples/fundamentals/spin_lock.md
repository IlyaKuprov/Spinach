# examples/fundamentals/spin_lock.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/spin_lock.m`
- Signature: `spin_lock()`
- Total lines: 66

## Purpose

A spin-locking experiment on a two-spin system.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- A spin-locking experiment on a two-spin system.
- Isotopes
- Magnetic induction
- Chemical shifts
- Scalar couplings
- Basis set
- Spinach housekeeping
- Initial state
- Observable states
- Pulse operators
- Hamiltonian
- Spin-locking field of 1.5 kHz along Y

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `step()`, `evolution()`, `kfigure()`, `plot3()`, `answer()`.
