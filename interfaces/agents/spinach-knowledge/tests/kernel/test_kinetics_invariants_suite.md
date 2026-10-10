# tests/kernel/test_kinetics_invariants_suite.m

## Purpose

Regression test for deterministic chemical kinetics helpers in Spinach. The suite verifies closed-form steady states, independent reaction blocks, first-order exchange routing and conservation, and empty reaction maps on a reordered local descriptor.

## Behaviour

The function announces the test target with `fprintf('TESTING: Chemical kinetics invariants\n')` and initialises a result object via `new_test_result('kernel/kinetics_invariants_suite', 'Chemical kinetics invariants', 'supported kinetic generators must conserve and route spin order.')`.

It then performs the following regression checks with explanatory messages:

- **Two-site steady state from detailed balance**: with `kf=2`, `kr=5`, `K=[-kf kr; kf -kr]`, `c0=[2;1]`, `ctot=sum(c0)`, and reference `c_ref=ctot*[kr; kf]/(kf+kr)`, it compares `equilibrate(K,c0)` to `c_ref` with tolerances `1e-13` (absolute and relative). The stated invariant: at equilibrium `k_forward c_1` equals `k_reverse c_2` and total concentration is conserved.
- **Independent reaction blocks**: with `K1=[-1 4;1 -4]`, `K2=[-3 2;3 -2]`, `K=blkdiag(K1,K2)`, `c0=[3;0;1;2]`, and reference `c_ref=[sum(c0(1:2))*[4;1]/5; sum(c0(3:4))*[2;3]/5]`, it compares `equilibrate(K,c0)` to `c_ref` with tolerances `1e-13`. The stated invariant: independent kinetic components equilibrate separately and retain their own material totals.
- **Zero-concentration shortcut**: compares `equilibrate(K,zeros(4,1))` to `zeros(4,1)` with tolerances `1e-15`. The stated invariant: a zero initial concentration vector remains zero for linear kinetics.
- **Exchange column sums**: a two-substance system transfers spin order from spin 1 to spin 2 at rate 3. The kinetic generator has zero column sums, checked with absolute and relative whole-vector tolerances of `1e-14`.
- **Flux routing**: the generator applied to `Lz` on spin 1 equals three times the difference between the destination and source states; its action on destination `Lz` is zero. Both comparisons use absolute and relative whole-vector tolerances of `1e-14`.
- **Reordered local descriptor**: with `chem.parts={[2 1]}`, an empty reaction returns an empty cell array. The empty-record check is independent of the preceding two-substance routing fixture.

## Inputs and outputs

- **Outputs**: `result` — regression test result with explanatory messages, returned by the function.
- The function takes no inputs.

## References

- Source: [tests/kernel/test_kinetics_invariants_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_kinetics_invariants_suite.m)

The one-way reaction record explicitly matches spin 1 to spin 2 and preserves the original routing and zero-destination-action assertions.
