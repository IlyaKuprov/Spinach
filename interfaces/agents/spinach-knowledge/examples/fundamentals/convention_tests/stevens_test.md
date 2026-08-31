# examples/fundamentals/convention_tests/stevens_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/convention_tests/stevens_test.m`
- Signature: `stevens_test()`
- Total lines: 123

## Purpose

Tests of Spinach Stevens operator function against explicit expressions from the literature.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Spin and multiplicity; implemented by `sqn=12; mult=2*sqn+1`.
- Lines 11-12: Spinach results, rank 6; implemented by `O_6_p6_Spinach=stevens(mult,6,+6)`.
- Lines 26-27: Explicit formulae, rank 6; implemented by `S=pauli(mult); s=sqn*(sqn+1)`.
- Lines 42-49: Differences, rank 6; implemented by `diffs=[norm(O_6_p6_Spinach-O_6_p6_Stevens,1) norm(O_6_p5_Spinach-O_6_p5_Stevens,1) norm(O_6_p4_Spinach-O_6_p4_Stevens,1) norm(O_6_p3_Spinach-O_6_p3_Stevens,1) norm(O_6_p…`.
- Lines 52-53: Spinach results, rank 4; implemented by `O_4_p4_Spinach=stevens(mult,4,+4)`.
- Lines 63-64: Explicit formulae, rank 4; implemented by `S=pauli(mult); s=sqn*(sqn+1)`.
- Lines 75-80: Differences, rank 4; implemented by `diffs=[norm(O_4_p4_Spinach-O_4_p4_Stevens,1) norm(O_4_p3_Spinach-O_4_p3_Stevens,1) norm(O_4_p2_Spinach-O_4_p2_Stevens,1) norm(O_4_p1_Spinach-O_4_p1_Stevens,1) norm(O_4_0…`.
- Lines 83-84: Spinach results, rank 2; implemented by `O_2_p2_Spinach=stevens(mult,2,+2)`.
- Lines 90-91: Explicit formulae, rank 2; implemented by `S=pauli(mult); s=sqn*(sqn+1)`.
- Lines 98-101: Differences, rank 2; implemented by `diffs=[norm(O_2_p2_Spinach-O_2_p2_Stevens,1) norm(O_2_p1_Spinach-O_2_p1_Stevens,1) norm(O_2_0_Spinach-O_2_0_Stevens,1) norm(O_2_m1_Spinach-O_2_m1_Stevens,1) norm(O_2_m2_…`.
- Lines 104-105: IST formulae, rank 2; implemented by `T=irr_sph_ten(mult,2)`.
- Lines 113-116: Differences, rank 2; implemented by `diffs=[norm(O_2_p2_Spinach-O_2_p2_IST,1) norm(O_2_p1_Spinach-O_2_p1_IST,1) norm(O_2_0_Spinach-O_2_0_IST,1) norm(O_2_m1_Spinach-O_2_m1_IST,1) norm(O_2_m2_Spinach-O_2_m2_I…`.
- Lines 119-120: Report that the test is passed; implemented by `disp('Stevens operator test PASSED.')`.

### Control flow inferred from the code

- Line 50: conditional branch on `norm(diffs)>1e-4, error('Stevens operator test FAILED.'); end`.
- Line 81: conditional branch on `norm(diffs)>1e-7, error('Stevens operator test FAILED.'); end`.
- Line 102: conditional branch on `norm(diffs)>1e-12, error('Stevens operator test FAILED.'); end`.
- Line 117: conditional branch on `norm(diffs)>1e-10, error('Stevens operator test FAILED.'); end`.

### Key state/data transformations

- Lines 9: computes `sqn` using `sqn=12; mult=2*sqn+1`.
- Lines 12: computes `O_6_p6_Spinach` using `O_6_p6_Spinach=stevens(mult,6,+6)`.
- Lines 13: computes `O_6_p5_Spinach` using `O_6_p5_Spinach=stevens(mult,6,+5)`.
- Lines 14: computes `O_6_p4_Spinach` using `O_6_p4_Spinach=stevens(mult,6,+4)`.
- Lines 15: computes `O_6_p3_Spinach` using `O_6_p3_Spinach=stevens(mult,6,+3)`.
- Lines 16: computes `O_6_p2_Spinach` using `O_6_p2_Spinach=stevens(mult,6,+2)`.
- Lines 17: computes `O_6_p1_Spinach` using `O_6_p1_Spinach=stevens(mult,6,+1)`.
- Lines 18: computes `O_6_0_Spinach` using `O_6_0_Spinach =stevens(mult,6, 0)`.
- Lines 19: computes `O_6_m1_Spinach` using `O_6_m1_Spinach=stevens(mult,6,-1)`.
- Lines 20: computes `O_6_m2_Spinach` using `O_6_m2_Spinach=stevens(mult,6,-2)`.
- Lines 21: computes `O_6_m3_Spinach` using `O_6_m3_Spinach=stevens(mult,6,-3)`.
- Lines 22: computes `O_6_m4_Spinach` using `O_6_m4_Spinach=stevens(mult,6,-4)`.
- Lines 23: computes `O_6_m5_Spinach` using `O_6_m5_Spinach=stevens(mult,6,-5)`.
- Lines 24: computes `O_6_m6_Spinach` using `O_6_m6_Spinach=stevens(mult,6,-6)`.
- Lines 27: computes `S` using `S=pauli(mult); s=sqn*(sqn+1)`.
- Lines 28: computes `O_6_p6_Stevens` using `O_6_p6_Stevens=(1/2)*(S.p^6+S.m^6)`.
- Lines 29: computes `O_6_p5_Stevens` using `O_6_p5_Stevens=(1/2)*acomm(S.z,S.p^5+S.m^5)/2`.
- Lines 30: computes `O_6_p4_Stevens` using `O_6_p4_Stevens=(1/2)*acomm(11*S.z^2-(s+38)*S.u,S.p^4+S.m^4)/2`.

## Implementation structure

- Tests of Spinach Stevens operator function against
- explicit expressions from the literature.
- Spin and multiplicity
- Spinach results, rank 6
- Explicit formulae, rank 6
- Differences, rank 6
- Spinach results, rank 4
- Explicit formulae, rank 4
- Differences, rank 4
- Spinach results, rank 2
- Explicit formulae, rank 2
- Differences, rank 2

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `stevens()`, `pauli()`, `acomm()`, `irr_sph_ten()`.
