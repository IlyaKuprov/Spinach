# examples/fundamentals/operator_tests/commutation_8.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/operator_tests/commutation_8.m`
- Signature: `commutation_8()`
- Total lines: 59

## Purpose

Delicate action and commutation tests for Hilbert-Liouville conversion and Stevens operators.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Accuracy threshold; implemented by `tol=1e-10`.
- Lines 11-12: Test Hilbert-Liouville conversion identities; implemented by `dim=6`.
- Lines 21-22: Build Hilbert-space products for vectorization checks; implemented by `HR=H*R; RH=R*H`.
- Lines 24-25: Compare left, right, commutator and anticommutator actions; implemented by `err_left=norm(L_left*R_vec-HR(:),2)`.
- Lines 30-31: Report Hilbert-Liouville conversion failures; implemented by `if max([err_left err_right err_comm err_acomm])>tol`.
- Lines 36-37: Test first-rank Stevens operators against Pauli operators; implemented by `for mult=[2 3 5 8]`.
- Lines 39-40: Build Stevens and Pauli operators; implemented by `L=pauli(mult)`.
- Lines 45-46: Compare first-rank operators; implemented by `err_z=norm(O_10-L.z,'fro')`.
- Lines 50-51: Report first-rank Stevens failures; implemented by `if max([err_z err_x err_y])>tol`.

### Control flow inferred from the code

- Line 31: conditional branch on `max([err_left err_right err_comm err_acomm])>tol`.
- Line 37: `for` loop over `mult=[2 3 5 8]`.
- Line 51: conditional branch on `max([err_z err_x err_y])>tol`.

### Key state/data transformations

- Lines 9: computes `tol` using `tol=1e-10`.
- Lines 12: computes `dim` using `dim=6`.
- Lines 13: computes `H` using `H=randn(dim)+1i*randn(dim)`.
- Lines 14: computes `R` using `R=randn(dim)+1i*randn(dim)`.
- Lines 15: computes `L_left` using `L_left=hilb2liouv(H,'left')`.
- Lines 16: computes `L_right` using `L_right=hilb2liouv(H,'right')`.
- Lines 17: computes `L_comm` using `L_comm=hilb2liouv(H,'comm')`.
- Lines 18: computes `L_acomm` using `L_acomm=hilb2liouv(H,'acomm')`.
- Lines 19: computes `R_vec` using `R_vec=hilb2liouv(R,'statevec')`.
- Lines 22: computes `HR` using `HR=H*R; RH=R*H`.
- Lines 25: computes `err_left` using `err_left=norm(L_left*R_vec-HR(:),2)`.
- Lines 26: computes `err_right` using `err_right=norm(L_right*R_vec-RH(:),2)`.
- Lines 27: computes `err_comm` using `err_comm=norm(L_comm*R_vec-reshape(HR-RH,[],1),2)`.
- Lines 28: computes `err_acomm` using `err_acomm=norm(L_acomm*R_vec-reshape(HR+RH,[],1),2)`.
- Lines 40: computes `L` using `L=pauli(mult)`.
- Lines 41: computes `O_10` using `O_10=stevens(mult,1,0)`.
- Lines 42: computes `O_11` using `O_11=stevens(mult,1,1)`.
- Lines 43: computes `O_1m1` using `O_1m1=stevens(mult,1,-1)`.

## Implementation structure

- Delicate action and commutation tests for Hilbert-Liouville
- conversion and Stevens operators.
- Accuracy threshold
- Test Hilbert-Liouville conversion identities
- Build Hilbert-space products for vectorization checks
- Compare left, right, commutator and anticommutator actions
- Report Hilbert-Liouville conversion failures
- Test first-rank Stevens operators against Pauli operators
- Build Stevens and Pauli operators
- Compare first-rank operators
- Report first-rank Stevens failures

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `hilb2liouv()`, `pauli()`, `stevens()`.
