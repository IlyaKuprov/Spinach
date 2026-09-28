# tests/kernel/test_dynamic_fp_contexts.m

- Signature: `result=test_dynamic_fp_contexts()`

## Purpose

Regression test for the `imaging()` and `meshflow()` context hand-off paths. Both contexts run a probe pulse sequence that returns their assembled operators and phantom-derived states for inspection.

## Physical / mathematical content

- The test checks generators acting on a combined spatial-and-spin state of size `spc_dim*spn_dim`.
- It checks conservation of total spatial mass through zero column sums of the flow generator: for periodic one-dimensional imaging flow and diffusion, and for finite-volume diffusion on a two-cell mesh with closed boundaries.

## Numerical / algorithmic content

- The imaging variant uses a one-spin `1H` spherical-tensor Liouville-space system on a ten-point, one-dimensional grid of extent `0.01`, with `deriv={'period',3}`, uniform velocity `1e-3`, and diffusion coefficient `1e-9`. It supplies relaxation, initial-state, and coil phantoms, then checks generator dimensions and finiteness, phantom-state lengths, and empty transverse gradient operators.
- The meshflow variant uses the same spin-system setup on a two-cell Voronoi mesh with a shared boundary, zero velocity, and diffusion coefficient `1e-8`. It supplies operator and state phantoms, then checks generator dimensions and finiteness and phantom-state lengths.
- The probe pulse sequence returns `H`, `R`, `K`, `G`, `F`, `rho0`, `coil`, and the spin, spatial, and combined dimensions. The matrix checks reject `NaN` and `Inf` values.

## Outputs

- `result` — regression test result with explanatory messages.

## Attribution

- ilya.kuprov@weizmann.ac.il