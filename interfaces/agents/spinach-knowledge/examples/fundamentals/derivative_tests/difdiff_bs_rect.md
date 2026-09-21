# examples/fundamentals/derivative_tests/difdiff_bs_rect.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/derivative_tests/difdiff_bs_rect.m`
- Signature: `difdiff_bs_rect()`
- Total lines: 151

## Purpose

Directional derivative test for Cartesian GRAPE with Bloch-Siegert corrections.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Implementation structure

- Directional derivative test for Cartesian GRAPE with Bloch-Siegert
- corrections.
- Set the magnetic field
- Set isotopes
- Set interactions
- Set basis
- Run Spinach housekeeping
- Build and normalise initial state
- Build and normalise target state
- Get control operators
- Get offset and shift operators
- Build drift Hamiltonian

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `singlet()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `true()`, `optimcon()`, `grape_xy()`, `squeeze()`, `grad_anl()`, `wave_forw()`, `wave_back()`, `fid_forw()`, `fid_back()`.
