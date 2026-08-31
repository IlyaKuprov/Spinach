# examples/fundamentals/state_tests/state_consistency_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_tests/state_consistency_1.m`
- Signature: `state_consistency_1()`
- Total lines: 76

## Purpose

Test of internal consistency for state and operator generation across the three formalisms supported by Spinach. Two-and four-spin states are tested.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 12-13: Set the spin system; implemented by `sys.isotopes={'1H','1H','1H','1H'}`.
- Lines 16-17: Loop over the formalisms; implemented by `Fs={'zeeman-hilb','zeeman-liouv','sphten-liouv'}`.
- Lines 20-21: Basis set; implemented by `bas.formalism=Fs{n}`.
- Lines 24-25: Hush the logs; implemented by `sys.output='hush'`.
- Lines 28-29: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Unit state from Spinach; implemented by `U=state(spin_system,{'E','E','E','E'},{1 2 3 4})`.
- Lines 35-37: Two-spin singlet-triplet state sum test; implemented by `combos={[1 2],[1 3],[1 4],[2 1],[2 3],[2 4], [3 1],[3 2],[3 4],[4 1],[4 2],[4 3]}`.
- Lines 46-47: Four-spin singlet-triplet state sum test; implemented by `combos=num2cell(perms([1 2 3 4]),2)`.
- Lines 70-71: Good news to the user; implemented by `disp([Fs{n} ': state construction test PASSED.'])`.

### Control flow inferred from the code

- Line 18: `for` loop over `n=1:numel(Fs)`.
- Line 38: `parfor` loop over `k=1:numel(combos)`.
- Line 41: conditional branch on `norm(S+TU+T0+TD-U,1)>1e-6`.
- Line 48: `parfor` loop over `k=1:numel(combos)`.
- Line 65: conditional branch on `norm(S-U,1)>1e-6`.

### Key state/data transformations

- Lines 10: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H'}`.
- Lines 14: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0 0 0 0}`.
- Lines 17: computes `Fs` using `Fs={'zeeman-hilb','zeeman-liouv','sphten-liouv'}`.
- Lines 21: computes `bas.formalism` using `bas.formalism=Fs{n}`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `sys.output` using `sys.output='hush'`.
- Lines 26: computes `sys.disable` using `sys.disable={'hygiene'}`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `U` using `U=state(spin_system,{'E','E','E','E'},{1 2 3 4})`.
- Lines 36-37: computes `combos` using `combos={[1 2],[1 3],[1 4],[2 1],[2 3],[2 4], [3 1],[3 2],[3 4],[4 1],[4 2],[4 3]}`.
- Lines 39: computes `S` using `S=singlet(spin_system,combos{k}(1),combos{k}(2))`.
- Lines 40: computes `[TU,T0,TD]` using `[TU,T0,TD]=triplet(spin_system,combos{k}(1),combos{k}(2))`.

## Implementation structure

- Test of internal consistency for state and operator
- generation across the three formalisms supported by
- Spinach. Two-and four-spin states are tested.
- Magnet field
- Set the spin system
- Loop over the formalisms
- Basis set
- Hush the logs
- Spinach housekeeping
- Unit state from Spinach
- Two-spin singlet-triplet state sum test
- Four-spin singlet-triplet state sum test

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `singlet()`, `triplet()`, `num2cell()`, `perms()`, `four_spin_states()`.
