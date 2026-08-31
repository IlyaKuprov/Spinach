# tests/kernel/test_operator_expansion_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_operator_expansion_suite.m`
- Signature: `result=test_operator_expansion_suite()`
- Total lines: 156

## Purpose

Tests operator expansion and conversion helpers. Syntax: result=test_operator_expansion_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file also defines local helper function(s): `ist_reconstruct()`, `bm_reconstruct()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Operator expansion and conversion helpers\n')`.
- Lines 20-23: State the expansion target of the test; implemented by `result=new_test_result('kernel/operator_expansion_suite', 'Operator expansions and conversions', 'operator expansion coefficients must reconstruct the source matrices.')`.
- Lines 25-26: Check Hilbert-to-Liouville vectorisation identities on a non-diagonal matrix; implemented by `H=[1 2;3 4]`.
- Lines 35-36: Check IST expansion of a generic spin-one matrix; implemented by `A=[1 2 0;0 -1 3;4 0 2]`.
- Lines 41-42: Check spin and boson energy-level counting conventions in IST expansions; implemented by `[states,coeffs]=enlev2ist(3,1,'S')`.
- Lines 51-52: Check central-transition and boson-product IST expansion wrappers; implemented by `[states,coeffs]=ct2ist(4,'+')`.
- Lines 60-61: Check bosonic monomial expansion of a generic finite oscillator matrix; implemented by `A=[1 2 0;0 -1 3;4 0 2]`.
- Lines 70-71: Build one-spin systems in the main formalisms; implemented by `sys.magnet=0`.
- Lines 76-77: Check Hilbert-space unit and sparse preallocation dimensions; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 85-86: Check Zeeman-Liouville unit and sparse preallocation dimensions; implemented by `bas.formalism='zeeman-liouv'`.
- Lines 94-95: Check spherical-tensor Liouville dimensions without assuming a fixed basis size literal; implemented by `bas.formalism='sphten-liouv'`.
- Lines 104-105: Check two-spin irreducible tensor formula in Hilbert space; implemented by `sys.isotopes={'1H','1H'}`.
- Lines 118-119: Check Lindbladian rate calibration on a diagonal jump process; implemented by `rho=[1;1]`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/operator_expansion_suite', 'Operator expansions and conversions', 'operator expansion coefficients must reconstruct the source matrices.')`.
- Lines 26: computes `H` using `H=[1 2;3 4]`.
- Lines 27: computes `unit` using `unit=speye(2)`.
- Lines 36: computes `A` using `A=[1 2 0;0 -1 3;4 0 2]`.
- Lines 37: computes `[states,coeffs]` using `[states,coeffs]=oper2ist(A)`.
- Lines 43: computes `P` using `P=zeros(3); P(3,3)=1`.
- Lines 55: computes `W` using `W=weyl(3)`.
- Lines 71: computes `sys.magnet` using `sys.magnet=0`.
- Lines 72: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 73: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 74: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 77: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 78: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 97: computes `basis_dim` using `basis_dim=size(spin_system.bas.basis,1)`.
- Lines 107: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=0`.
- Lines 108: computes `inter.coupling.scalar{2,2}` using `inter.coupling.scalar{2,2}=0`.
- Lines 111: computes `T20` using `T20=twospinist(spin_system,1,2,[2 0],'comm')`.
- Lines 112-114: computes `T20_ref` using `T20_ref=sqrt(2/3)*(operator(spin_system,{'Lz','Lz'},{1,2})- (1/4)*(operator(spin_system,{'L+','L-'},{1,2})+ operator(spin_system,{'L-','L+'},{1,2})))`.

### Local helper functions

- Line 129: `ist_reconstruct()` — `function A=ist_reconstruct(mult,states,coeffs)`. Build the complete IST basis for the requested multiplicity
  - Representative operation: `T=irr_sph_ten(mult)`.
  - Representative operation: `A=zeros(mult)`.
- Line 143: `bm_reconstruct()` — `function A=bm_reconstruct(nlevels,states,coeffs)`. Build the complete bosonic monomial basis for the requested truncation
  - Representative operation: `B=boson_mono(nlevels)`.
  - Representative operation: `A=zeros(nlevels)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks that irreducible spherical tensor and bosonic monomial
- expansion helpers reconstruct explicit matrices, and that operator-sized
- allocation helpers return the correct formalism dimensions.

## Implementation structure

- Tests operator expansion and conversion helpers. Syntax:
- result=test_operator_expansion_suite()
- result -regression test result with explanatory messages
- The test checks that irreducible spherical tensor and bosonic monomial
- expansion helpers reconstruct explicit matrices, and that operator-sized
- allocation helpers return the correct formalism dimensions.
- Announce the test target
- State the expansion target of the test
- Check Hilbert-to-Liouville vectorisation identities on a non-diagonal matrix
- Check IST expansion of a generic spin-one matrix
- Check spin and boson energy-level counting conventions in IST expansions
- Check central-transition and boson-product IST expansion wrappers

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `speye()`, `test_close()`, `hilb2liouv()`, `transpose()`, `oper2ist()`, `ist_reconstruct()`, `enlev2ist()`, `ct2ist()`, `centrans()`, `centran()`, `weyl()`, `bos2ist()`, `oper2bm()`, `bm_reconstruct()`, `enlev2bm()`.
