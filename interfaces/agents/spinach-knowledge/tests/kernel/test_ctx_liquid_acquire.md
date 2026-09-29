# tests/kernel/test_ctx_liquid_acquire.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_ctx_liquid_acquire.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_ctx_liquid_acquire.m)

## Purpose

Regression test that verifies the liquid context (`liquid()`) against the direct `acquire()` path. The test checks that `liquid()` passes the offset Liouvillian to `acquire()` correctly.

## Behaviour

1. Announces the test target with `fprintf('TESTING: Liquid context acquire path\n')` and initialises a test result via `new_test_result` for `kernel/ctx_liquid_acquire`.
2. Builds a one-spin Liouville-space spin system:
   - `sys.magnet = 14.1` (field in Tesla)
   - `sys.isotopes = {'1H'}`
   - `inter.zeeman.scalar = {0}`
   - `bas.formalism = 'sphten-liouv'`
   - `bas.approximation = 'none'`
   - `bas.projections = {+1}`
3. Sets up a short offset acquisition:
   - `parameters.spins = {'1H'}`
   - `parameters.rho0 = state(spin_system,'L+','1H')`
   - `parameters.coil = state(spin_system,'L+','1H')`
   - `parameters.decouple = {}`
   - `parameters.offset = 25`
   - `parameters.sweep = 1000`
   - `parameters.npoints = 4`
   - `parameters.needs = {}`
   - `parameters.rframes = {}`
4. Runs the production liquid context: `fid_ctx = liquid(spin_system,@acquire,parameters,'nmr')`.
5. Builds the reference direct path: applies `assume(spin_system,'nmr')`, computes `H = hamiltonian(spin_system)`, `R = relaxation(spin_system)`, `K = kinetics(spin_system)`, applies `H = frqoffset(spin_system,H,parameters)`, then calls `fid_ref = acquire(spin_system,parameters,H,R,K)`.
6. Compares the two FIDs with `test_close(result,'liquid context FID',fid_ctx,fid_ref,1e-12,1e-12,...)`, asserting that `liquid()` reproduces direct `acquire()` for the same offset Liouvillian.

## Inputs and outputs

```matlab
result = test_ctx_liquid_acquire()
```

- **Output:** `result` — regression test result structure with explanatory messages.
- **Input:** none.

## References
