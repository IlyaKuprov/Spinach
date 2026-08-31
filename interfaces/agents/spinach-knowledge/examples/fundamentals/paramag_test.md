# examples/fundamentals/paramag_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/paramag_test.m`
- Signature: `paramag_test()`
- Total lines: 32

## Purpose

Paramagnetic chemical shift module tests.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 7-8: xyz2pms against ppcs; implemented by `chi=randn(3); chi=(chi+chi')/2`.
- Lines 18-19: xyz2hfc + hfc2pms against xyz2pms; implemented by `chi=randn(3); chi=(chi+chi')/2`.

### Control flow inferred from the code

- Line 12: conditional branch on `abs(pcs_a-pcs_b)<1e-6`.
- Line 25: conditional branch on `norm(pms_tensor_a-pms_tensor_b,'fro')<1e-6`.

### Key state/data transformations

- Lines 8: computes `chi` using `chi=randn(3); chi=(chi+chi')/2`.
- Lines 9: computes `nxyz` using `nxyz=randn(1,3); sxyz=randn(1,3)`.
- Lines 10: computes `pcs_a` using `pcs_a=ppcs(nxyz,sxyz,chi)`.
- Lines 11: computes `pcs_b` using `pcs_b=trace(xyz2pms(nxyz,sxyz,chi))/3`.
- Lines 21: computes `isotope` using `isotope='15N'`.
- Lines 22: computes `A` using `A=xyz2hfc(mxyz,nxyz,isotope)`.
- Lines 23: computes `[~,pms_tensor_a]` using `[~,pms_tensor_a]=hfc2pms(A,chi,isotope)`.
- Lines 24: computes `pms_tensor_b` using `pms_tensor_b=xyz2pms(nxyz,mxyz,chi)`.

## Implementation structure

- Paramagnetic chemical shift module tests.
- xyz2pms against ppcs
- xyz2hfc + hfc2pms against xyz2pms

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `ppcs()`, `xyz2pms()`, `xyz2hfc()`, `hfc2pms()`.
