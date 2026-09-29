# tests/kernel/test_kinetics_invariants_suite.m

## Purpose

Regression test for deterministic chemical kinetics helpers in Spinach. The suite verifies closed-form steady states, independent reaction blocks, reaction-generator state routing, and conservation in a tiny exchange kinetics superoperator.

## Behaviour

The function announces the test target with `fprintf('TESTING: Chemical kinetics invariants\n')` and initialises a result object via `new_test_result('kernel/kinetics_invariants_suite', 'Chemical kinetics invariants', 'kinetic generators must conserve matter and route spin order between declared species.')`.

It then performs the following checks, each through `test_close` with explanatory messages:

- **Two-site steady state from detailed balance**: with `kf=2`, `kr=5`, `K=[-kf kr; kf -kr]`, `c0=[2;1]`, `ctot=sum(c0)`, and reference `c_ref=ctot*[kr; kf]/(kf+kr)`, it compares `equilibrate(K,c0)` to `c_ref` with tolerances `1e-13` (absolute and relative). The stated invariant: at equilibrium `k_forward c_1` equals `k_reverse c_2` and total concentration is conserved.
- **Independent reaction blocks**: with `K1=[-1 4;1 -4]`, `K2=[-3 2;3 -2]`, `K=blkdiag(K1,K2)`, `c0=[3;0;1;2]`, and reference `c_ref=[sum(c0(1:2))*[4;1]/5; sum(c0(3:4))*[2;3]/5]`, it compares `equilibrate(K,c0)` to `c_ref` with tolerances `1e-13`. The stated invariant: independent kinetic components equilibrate separately and retain their own material totals.
- **Zero-concentration shortcut**: compares `equilibrate(K,zeros(4,1))` to `zeros(4,1)` with tolerances `1e-15`. The stated invariant: a zero initial concentration vector remains zero for linear kinetics.
- **Closed-exchange column sums**: builds a two-site spherical-tensor spin system with `sys.magnet=14.1`, `sys.isotopes={'1H','1H'}`, `inter.zeeman.scalar={0,0}`, `inter.chem.parts={1,2}`, `inter.chem.rates=[-3 3;3 -3]`, `inter.chem.concs=[1 1]`, `bas.formalism='sphten-liouv'`, `bas.approximation='none'`, then `spin_system=test_spin_system(sys,inter,bas)`. It computes `K=kinetics(spin_system)` and `col_sums=sum(full(K),1)`, comparing to zeros with tolerances `1e-14`. The stated invariant: in a closed two-site exchange model all probability leaving a column re-enters elsewhere.
- **Reaction-generator routing**: with `reaction.reactants=1`, `reaction.products=2`, `reaction.matching=[1 2]`, it computes `G=react_gen(spin_system,reaction)`, `rho_source=state(spin_system,'Lz',1)`, `rho_destin=state(spin_system,'Lz',2)`. It checks `G{1}*rho_source` equals `rho_destin-rho_source` (tolerances `1e-14`), i.e. a matched reactant spin order is removed from the source species and inserted into the product species; and `G{1}*rho_destin` equals zeros (tolerances `1e-14`), i.e. the reactant generator does not drain states already located on the product species.

## Inputs and outputs

- **Outputs**: `result` — regression test result with explanatory messages, returned by the function.
- The function takes no inputs.

## References

- Source: [tests/kernel/test_kinetics_invariants_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_kinetics_invariants_suite.m)
