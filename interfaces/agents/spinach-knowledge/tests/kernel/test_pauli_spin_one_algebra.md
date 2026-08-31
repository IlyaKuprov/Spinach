# tests/kernel/test_pauli_spin_one_algebra.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_pauli_spin_one_algebra.m`
- Signature: `result=test_pauli_spin_one_algebra()`
- Total lines: 52

## Purpose

Tests spin-one angular momentum matrices. Syntax: result=test_pauli_spin_one_algebra()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Spin-one angular momentum algebra\n')`.
- Lines 19-22: State the physical target of the test; implemented by `result=new_test_result('kernel/pauli_spin_one_algebra', 'Spin-one angular momentum algebra', 'Spin-one operators must realise su(2), with S^2=2.')`.
- Lines 24-25: Generate Spinach spin-one operators; implemented by `S=pauli(3)`.
- Lines 27-28: Write the textbook spin-one matrices explicitly; implemented by `Sp=[0 sqrt(2) 0;0 0 sqrt(2);0 0 0]`.
- Lines 34-36: Check matrix elements and commutators; implemented by `result=test_close(result,'Sz projections',S.z,Sz,1e-15,1e-15, 'spin-one magnetic quantum numbers are +1, 0, and -1')`.
- Lines 46-47: Check the Casimir operator; implemented by `S2=S.x*S.x+S.y*S.y+S.z*S.z`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/pauli_spin_one_algebra', 'Spin-one angular momentum algebra', 'Spin-one operators must realise su(2), with S^2=2.')`.
- Lines 22: computes `'Spin-one operators must realise su(2), with S^2` using `'Spin-one operators must realise su(2), with S^2=2.')`.
- Lines 25: computes `S` using `S=pauli(3)`.
- Lines 28: computes `Sp` using `Sp=[0 sqrt(2) 0;0 0 sqrt(2);0 0 0]`.
- Lines 29: computes `Sm` using `Sm=Sp'`.
- Lines 30: computes `Sz` using `Sz=diag([1 0 -1])`.
- Lines 31: computes `Sx` using `Sx=(Sp+Sm)/2`.
- Lines 32: computes `Sy` using `Sy=(Sp-Sm)/2i`.
- Lines 47: computes `S2` using `S2=S.x*S.x+S.y*S.y+S.z*S.z`.
- Lines 49: computes `'for s` using `'for s=1 the Casimir eigenvalue is 2')`.

## Outputs

- result -regression test result with explanatory messages
- The test checks the spin-one representation: Sz projections are +1, 0,
- and -1, ladder matrix elements are sqrt(2), and S^2=s(s+1)=2.

## Implementation structure

- Tests spin-one angular momentum matrices. Syntax:
- result=test_pauli_spin_one_algebra()
- result -regression test result with explanatory messages
- The test checks the spin-one representation: Sz projections are +1, 0,
- and -1, ladder matrix elements are sqrt(2), and S^2=s(s+1)=2.
- Announce the test target
- State the physical target of the test
- Generate Spinach spin-one operators
- Write the textbook spin-one matrices explicitly
- Check matrix elements and commutators
- Check the Casimir operator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `pauli()`, `test_close()`, `comm()`.
