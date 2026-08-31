# tests/kernel/test_chemical_exchange_conservation.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_chemical_exchange_conservation.m`
- Signature: `result=test_chemical_exchange_conservation()`
- Total lines: 44

## Purpose

Tests conservation in two-site chemical exchange. Syntax: result=test_chemical_exchange_conservation()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Chemical exchange conservation\n')`.
- Lines 19-22: State the kinetics target of the test; implemented by `result=new_test_result('kernel/chemical_exchange_conservation', 'Chemical exchange conservation', 'closed two-site exchange must conserve total spin population.')`.
- Lines 24-25: Build a symmetric two-site exchange system; implemented by `sys.magnet=14.1`.
- Lines 35-36: Build the kinetics generator; implemented by `K=kinetics(spin_system)`.
- Lines 38-39: Closed Markov kinetics conserve total population by zero column sums; implemented by `col_sums=sum(full(K),1)`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/chemical_exchange_conservation', 'Chemical exchange conservation', 'closed two-site exchange must conserve total spin population.')`.
- Lines 25: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0 0}`.
- Lines 28: computes `inter.chem.parts` using `inter.chem.parts={1,2}`.
- Lines 29: computes `inter.chem.rates` using `inter.chem.rates=[-3 3;3 -3]`.
- Lines 30: computes `inter.chem.concs` using `inter.chem.concs=[1 1]`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 33: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 36: computes `K` using `K=kinetics(spin_system)`.
- Lines 39: computes `col_sums` using `col_sums=sum(full(K),1)`.

## Outputs

- result -regression test result with explanatory messages
- The test builds a symmetric two-site exchange model and checks that the
- kinetics generator conserves the total population over the two sites.

## Implementation structure

- Tests conservation in two-site chemical exchange. Syntax:
- result=test_chemical_exchange_conservation()
- result -regression test result with explanatory messages
- The test builds a symmetric two-site exchange model and checks that the
- kinetics generator conserves the total population over the two sites.
- Announce the test target
- State the kinetics target of the test
- Build a symmetric two-site exchange system
- Build the kinetics generator
- Closed Markov kinetics conserve total population by zero column sums

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `kinetics()`, `test_close()`.
