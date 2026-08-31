# tests/kernel/test_orientation_average_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_orientation_average_suite.m`
- Signature: `result=test_orientation_average_suite()`
- Total lines: 68

## Purpose

Tests orientation() and average() on small exact cases. Syntax: result=test_orientation_average_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Orientation and average helpers\n')`.
- Lines 20-23: State the rotational-kernel target of the test; implemented by `result=new_test_result('kernel/orientation_average_suite', 'Orientation and average helpers', 'rotational helper kernels must preserve exact limiting cases.')`.
- Lines 25-26: Build a synthetic rank-one rotational basis; implemented by `Q=cell(1,1)`.
- Lines 37-38: Contract the zero-orientation Hamiltonian; implemented by `H=orientation(Q,[0 0 0])`.
- Lines 41-43: Check the zero-angle Wigner identity path; implemented by `result=test_close(result,'zero Euler orientation',H,H_ref,1e-14,1e-14, 'wigner(r,0,0,0) is the identity, so only diagonal rotational components remain')`.
- Lines 45-46: Build a quiet spin system for average() diagnostics; implemented by `sys.magnet=0`.
- Lines 53-54: Define an unmodulated Hamiltonian decomposition; implemented by `Hp=sparse(2,2)`.
- Lines 59-60: Run the production average Hamiltonian helper; implemented by `H_avg=average(spin_system,Hp,H0,Hm,omega,'ah_first_order')`.
- Lines 62-64: Check the exact unmodulated limit; implemented by `result=test_close(result,'unmodulated average',H_avg,H0,1e-14,1e-14, 'with zero positive and negative Fourier components first-order averaging returns H0')`.

### Control flow inferred from the code

- Line 28: `for` loop over `n=1:3`.
- Line 29: `for` loop over `k=1:3`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/orientation_average_suite', 'Orientation and average helpers', 'rotational helper kernels must preserve exact limiting cases.')`.
- Lines 26: computes `Q` using `Q=cell(1,1)`.
- Lines 27: computes `Q{1}` using `Q{1}=cell(3,3)`.
- Lines 30: computes `Q{1}{n,k}` using `Q{1}{n,k}=sparse(2,2)`.
- Lines 33: computes `Q{1}{1,1}` using `Q{1}{1,1}=sparse([1 0;0 0])`.
- Lines 34: computes `Q{1}{2,2}` using `Q{1}{2,2}=sparse([0 1;1 0])`.
- Lines 35: computes `Q{1}{3,3}` using `Q{1}{3,3}=sparse([0 0;0 2])`.
- Lines 38: computes `H` using `H=orientation(Q,[0 0 0])`.
- Lines 39: computes `H_ref` using `H_ref=Q{1}{1,1}+Q{1}{2,2}+Q{1}{3,3}`.
- Lines 46: computes `sys.magnet` using `sys.magnet=0`.
- Lines 47: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 48: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 49: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 50: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 51: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 54: computes `Hp` using `Hp=sparse(2,2)`.
- Lines 55: computes `H0` using `H0=sparse([0 1;-1 0])`.
- Lines 56: computes `Hm` using `Hm=sparse(2,2)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks the zero-Euler-angle orientation contraction against an
- explicit diagonal Wigner sum and verifies that first-order average
- Hamiltonian theory leaves an unmodulated Hamiltonian unchanged.

## Implementation structure

- Tests orientation() and average() on small exact cases. Syntax:
- result=test_orientation_average_suite()
- result -regression test result with explanatory messages
- The test checks the zero-Euler-angle orientation contraction against an
- explicit diagonal Wigner sum and verifies that first-order average
- Hamiltonian theory leaves an unmodulated Hamiltonian unchanged.
- Announce the test target
- State the rotational-kernel target of the test
- Build a synthetic rank-one rotational basis
- Contract the zero-orientation Hamiltonian
- Check the zero-angle Wigner identity path
- Build a quiet spin system for average() diagnostics

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `orientation()`, `test_close()`, `wigner()`, `test_spin_system()`, `average()`.
