# tests/kernel/test_dynamic_fp_contexts.m

**Source**: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_fp_contexts.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_fp_contexts.m)

## Purpose

Regression test for the compact `imaging()` and `meshflow()` context hand-off paths. It verifies that both contexts assemble finite, correctly sized generators and phantom-derived initial and detection states, and that spatial flow is conserved.

## Behavior

- Announces the test target with `fprintf` and initializes a test result via `new_test_result` under the identifier `kernel/dynamic_fp_contexts`, describing the requirement that `imaging()` and `meshflow()` assemble context generators and phantom states with conserved spatial flow.
- Runs two subtests: `local_test_imaging` (Cartesian-grid imaging context) and `local_test_meshflow` (unstructured-mesh flow context).

### Imaging subtest

- Builds a one-spin `1H` spherical-tensor Liouville-space spin system (`bas.formalism='sphten-liouv'`, `bas.approximation='none'`, `sys.magnet=14.1`, zero scalar Zeeman interaction) via `test_spin_system`; the spin dimension is taken from `size(spin_system.bas.basis,1)`.
- Sets a minimal one-dimensional imaging grid: `parameters.npts=10`, `parameters.dims=0.01`, `parameters.offset=0`, `parameters.deriv={'period',3}`, `parameters.u=1e-3*ones(parameters.npts,1)`, `parameters.diff=1e-9`, with `parameters.spins={'1H'}` and `parameters.decouple={}`.
- Supplies relaxation, initial-state, and coil phantoms: `rlx_ph` zeros, `rlx_op` sparse `spn_dim`-by-`spn_dim`, `rho0_ph` ones, `rho0_st=state(spin_system,'Lz','1H')`, `coil_ph` ones, `coil_st=state(spin_system,'L+','1H')`.
- Calls `imaging(spin_system,@local_context_probe,parameters)`; the problem dimension is `parameters.npts*spn_dim`.
- Checks that `answer.spn_dim` equals the basis dimension and `answer.problem_dim` equals the Fokker-Planck size (Cartesian spatial dimension times spin dimension), both with zero tolerances.
- Checks via `local_test_square` that `answer.H`, `answer.R`, `answer.K`, `answer.F` (labelled `imaging Fx size`), and `answer.G{1}` (labelled `imaging Gx size`) are square matrices of the problem dimension and contain only finite values.
- Checks that `answer.G{2}` and `answer.G{3}` are both empty for a one-dimensional grid.
- Checks that `numel(answer.rho0)` and `numel(answer.coil)` equal the problem dimension (one spin/detection state per imaging voxel).
- Checks flow conservation: `sum(full(answer.F),1)` equals a zero row vector of length `problem_dim` with absolute and relative tolerances `1e-10`, reflecting that periodic one-dimensional flow and diffusion conserve total spatial mass.

### Meshflow subtest

- Builds the same one-spin `1H` spherical-tensor Liouville-space spin system; the spin dimension is `size(spin_system.bas.basis,1)`.
- Attaches a two-cell Voronoi mesh via `local_two_cell_mesh`: `mesh.idx.active=[1 2]`, `mesh.idx.triangles=[1 2 3]`, `mesh.vor.ncells=2`, `mesh.vor.weights=[1;1]`, `mesh.vor.cells={[1 2 3],[2 1 4]}`, `mesh.vor.vertices=[0 0;0 1;-1 0;1 0]`, cell-centre coordinates `mesh.x=[0;1]`, `mesh.y=[0;0]`, and zero velocity field `mesh.u=[0;0]`, `mesh.v=[0;0]`. The spatial dimension is `spin_system.mesh.vor.ncells` and the problem dimension is `spc_dim*spn_dim`.
- Supplies operator and state phantoms: `H_ph` ones, `H_op` sparse `spn_dim`-by-`spn_dim`, `R_ph` zeros, `R_op` sparse, `K_ph` zeros, `K_op=speye(spn_dim)`, `rho0_ph` ones, `rho0_st=state(spin_system,'Lz','1H')`, `coil_ph` ones, `coil_st=state(spin_system,'Lz','1H')`, and `parameters.diff=1e-8`.
- Calls `meshflow(spin_system,@local_context_probe,parameters)`.
- Checks that `answer.spn_dim` equals the basis dimension and `answer.problem_dim` equals the mesh cell count times spin dimension, both with zero tolerances.
- Checks via `local_test_square` that `answer.H`, `answer.R`, `answer.K`, and `answer.F` are square matrices of the problem dimension and contain only finite values.
- Checks that `numel(answer.rho0)` and `numel(answer.coil)` equal the problem dimension (one spin/detection state per active Voronoi cell).
- Checks flow conservation: `sum(full(answer.F),1)` equals a zero row vector of length `problem_dim` with absolute and relative tolerances `1e-15`, reflecting that finite-volume mesh diffusion with closed boundaries conserverves total mass.

### Helper functions

- `local_context_probe(~,parameters,H,R,K,G,F)` returns the context products for regression checks: `spc_dim`, `spn_dim`, `problem_dim` (as `parameters.spc_dim*parameters.spn_dim`), `rho0`, `coil`, `H`, `R`, `K`, `G`, and `F`.
- `local_test_square(result,label,A,matrix_dim)` uses `test_close` to verify `size(A)` equals `[matrix_dim matrix_dim]` with zero tolerances, and `test_true` to verify all elements of `full(A(:))` are finite (no NaN or Inf).

## Inputs and outputs

```matlab
result = test_dynamic_fp_contexts()
```

- **Output**: `result` — regression test result structure with explanatory messages, accumulated through `test_close` and `test_true` assertions.
- Takes no inputs.

## References

- `imaging` — Cartesian-grid imaging context exercised by this test.
- `meshflow` — unstructured-mesh flow context exercised by this test.
- `new_test_result`, `test_close`, `test_true`, `test_spin_system`, `state` — test harness and spin-system utilities used by this test.
