# examples/fundamentals/state_tests/state_consistency_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_tests/state_consistency_2.m`
- Signature: `state_consistency_2()`
- Total lines: 72

## Purpose

Test of consistency in the projection between spherical tensor basis set and Zeeman basis set.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Magneti field; implemented by `sys.magnet=14.1`.
- Lines 11-12: Isotopes; implemented by `sys.isotopes={'14N','235U'}`.
- Lines 14-15: No interactions; implemented by `inter.zeeman.scalar={0 0}`.
- Lines 17-18: Hilbert space, Zeeman basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 23-25: A suitably complicated state; implemented by `Op{1}=state(spin_system,{'Lz','Lx'},{1,2})+ state(spin_system,{'L+'},{1})`.
- Lines 27-28: Liouville space, Zeeman basis; implemented by `bas.formalism='zeeman-liouv'`.
- Lines 33-35: A suitably complicated state; implemented by `Op{2}=state(spin_system,{'Lz','Lx'},{1,2})+ state(spin_system,{'L+'},{1})`.
- Lines 37-38: Fold back into Hilbert space; implemented by `Op{2}=reshape(Op{2},[24 24])`.
- Lines 40-41: Liouville space, IST basis; implemented by `bas.formalism='sphten-liouv'`.
- Lines 46-48: A suitably complicated state; implemented by `Op{3}=state(spin_system,{'Lz','Lx'},{1,2})+ state(spin_system,{'L+'},{1})`.
- Lines 50-51: Project into Zeeman basis; implemented by `Op{3}=sphten2zeeman(spin_system)*Op{3}`.
- Lines 53-54: Fold back into Hilbert space; implemented by `Op{3}=reshape(Op{3},[24 24])`.
- Lines 56-59: Check the results; implemented by `if (norm(Op{1}-Op{2},1)>1e-6)|| (norm(Op{2}-Op{3},1)>1e-6)|| (norm(Op{3}-Op{1},1)>1e-6)`.
- Lines 61-62: Complain and bomb out; implemented by `error('State consistency test FAILED.')`.
- Lines 66-67: Good news to the user; implemented by `disp('State consistency test PASSED.')`.

### Control flow inferred from the code

- Line 57: conditional branch on `(norm(Op{1}-Op{2},1)>1e-6)||`.

### Key state/data transformations

- Lines 9: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'14N','235U'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0 0}`.
- Lines 18: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 19: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 20: computes `spin_system` using `spin_system=create(sys, inter)`.
- Lines 24-25: computes `Op{1}` using `Op{1}=state(spin_system,{'Lz','Lx'},{1,2})+ state(spin_system,{'L+'},{1})`.
- Lines 34-35: computes `Op{2}` using `Op{2}=state(spin_system,{'Lz','Lx'},{1,2})+ state(spin_system,{'L+'},{1})`.
- Lines 47-48: computes `Op{3}` using `Op{3}=state(spin_system,{'Lz','Lx'},{1,2})+ state(spin_system,{'L+'},{1})`.

## Implementation structure

- Test of consistency in the projection between spherical
- tensor basis set and Zeeman basis set.
- Magneti field
- Isotopes
- No interactions
- Hilbert space, Zeeman basis
- A suitably complicated state
- Liouville space, Zeeman basis
- Fold back into Hilbert space
- Liouville space, IST basis
- Project into Zeeman basis
- Check the results

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `sphten2zeeman()`.
