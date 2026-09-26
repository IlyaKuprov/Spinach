# examples/fundamentals/derivative_tests/dirdiff_test_system.m

- Signature: `[spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalism)`

## Purpose

Spin system generator for directional derivative tests. Syntax: [spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalism)

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Spin system generator for directional derivative tests. Syntax:
- [spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalism)
- Check consistency
- Select system size
- Keep the original large Liouville-space test
- Use a compact system for full Zeeman formalisms
- Set the magnetic field
- Put non-interacting spins at equal intervals
- within the [-100,+100] ppm chemical shift range
- Select the requested basis set
- Keep complete single-spin terms only
- Keep the full Zeeman basis
