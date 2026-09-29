# tests/kernel/test_kinetics_generator_suite.m

## Purpose

Regression test suite for the kinetics and flow generator helpers in Spinach. It verifies that `equilibrate`, `react_gen`, `kinetics`, and `flow_gen` satisfy conservation and detailed-balance invariants.

## Behaviour

- Announces the test target with `fprintf` and initialises a regression test result via `new_test_result` under the identifier `kernel/kinetics_generator_suite`, with the requirement that kinetic generators conserve matter and equilibrate closed systems correctly.
- **`equilibrate` two-state detailed balance**: builds the reversible Markov generator `K=[-2 1;2 -1]` with initial concentrations `c0=[3;0]`, and checks that the equilibrium concentrations equal `[1;2]` to absolute and relative tolerances `1e-14`, i.e. at equilibrium `k_21*c_1=k_12*c_2` while total concentration is conserved.
- **`equilibrate` zero shortcut**: checks that zero initial concentration `[0;0]` remains `[0;0]` with zero tolerances.
- Builds a two-site exchange spin system: `sys.magnet=0`, `sys.isotopes={'1H','1H'}`, `inter.zeeman.scalar={0 0}`, `inter.chem.parts={1,2}`, `inter.chem.rates=[-1 1;1 -1]`, `inter.chem.concs=[1 1]`, and `bas.formalism='sphten-liouv'` with `bas.approximation='none'`, passed through `test_spin_system`.
- **`react_gen` one reaction count**: with `reaction.reactants=1`, `reaction.products=2`, and `reaction.matching=[1 2]`, checks that `numel(G)` equals 1 with zero tolerances, since one reactant channel produces one reaction-generator matrix.
- **`react_gen` conservation**: checks that column sums of `full(G{1})` are zero (to `1e-14`), i.e. reaction drain and fill terms conserve total population column by column.
- **`kinetics` column sums**: checks that the full chemical kinetics generator `Kspin` for the symmetric exchange system has zero column sums (to `1e-14`), as expected for a closed chemical kinetics generator.
- **`flow_gen` minimal diffusion**: constructs a minimal two-cell mesh with `mesh.vor.ncells=2`, `mesh.vor.weights=[1;1]`, `mesh.vor.vertices=[0 0; 0 1]`, `mesh.vor.cells={[1 2],[1 2]}`, `mesh.idx.active=[1;2]`, `mesh.idx.triangles=[1 2 3]`, `mesh.x=[0;1;0]`, `mesh.y=[0;0;1]`, and zero velocity components `mesh.u`, `mesh.v`. With `flow_system.sys.output='hush'` and diffusion option `struct('diff',0.5)`, checks that `F` equals `[-0.5 0.5;0.5 -0.5]` to `1e-14`, since two identical cells sharing a unit boundary with `D=0.5` give a symmetric conservative diffusion generator.
- **`flow_gen` column sums**: checks that column sums of `full(F)` equal `[0 0]` to `1e-14`, i.e. the hydrodynamic flow/diffusion generator conserves total population.

## Inputs and outputs

- **Syntax**: `result=test_kinetics_generator_suite()`
- **Outputs**: `result` — regression test result with explanatory messages, accumulated through repeated `test_close` calls.
- **Inputs**: none.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_kinetics_generator_suite.m)
