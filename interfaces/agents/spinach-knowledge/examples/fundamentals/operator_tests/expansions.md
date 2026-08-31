# examples/fundamentals/operator_tests/expansions.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/operator_tests/expansions.m`
- Signature: `expansions()`
- Total lines: 83

## Purpose

Product action tests for irreducible spherical tensor operators and orthogonalised bosonic monomials.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Random level counts; implemented by `spin_mult=1+randi(9)`.
- Lines 12-13: Monomials and ISTs; implemented by `T=irr_sph_ten(spin_mult)`.
- Lines 16-17: Multiplication tables; implemented by `[ist_PTL,ist_PTR]=ist_product_table(spin_mult)`.
- Lines 20-21: IST test loop; implemented by `for n=1:numel(T)`.
- Lines 24-25: Left side product action; implemented by `LSPA_table=zeros(spin_mult,spin_mult)`.
- Lines 31-32: Right side product action; implemented by `RSPA_table=zeros(spin_mult,spin_mult)`.
- Lines 38-39: Accuracy tests; implemented by `norm_diff_a=norm(LSPA_table-LSPA_known,'fro')`.
- Lines 51-52: BM test loop; implemented by `for n=1:numel(B)`.
- Lines 55-56: Left side product action; implemented by `LSPA_table=zeros(bos_nlevels,bos_nlevels)`.
- Lines 62-63: Right side product action; implemented by `RSPA_table=zeros(bos_nlevels,bos_nlevels)`.

### Control flow inferred from the code

- Line 21: `for` loop over `n=1:numel(T)`.
- Line 22: `for` loop over `m=1:numel(T)`.
- Line 26: `for` loop over `k=1:numel(T)`.
- Line 33: `for` loop over `k=1:numel(T)`.
- Line 42: conditional branch on `(norm_diff_a/norm_scale>tol)||`.
- Line 52: `for` loop over `n=1:numel(B)`.
- Line 53: `for` loop over `m=1:numel(B)`.
- Line 57: `for` loop over `k=1:numel(B)`.
- Line 64: `for` loop over `k=1:numel(B)`.
- Line 73: conditional branch on `(norm_diff_a/norm_scale>tol)||`.

### Key state/data transformations

- Lines 9: computes `spin_mult` using `spin_mult=1+randi(9)`.
- Lines 10: computes `bos_nlevels` using `bos_nlevels=2+randi(8)`.
- Lines 13: computes `T` using `T=irr_sph_ten(spin_mult)`.
- Lines 14: computes `B` using `B=boson_ortho(bos_nlevels)`.
- Lines 17: computes `[ist_PTL,ist_PTR]` using `[ist_PTL,ist_PTR]=ist_product_table(spin_mult)`.
- Lines 18: computes `[bos_PTL,bos_PTR]` using `[bos_PTL,bos_PTR]=bos_product_table(bos_nlevels)`.
- Lines 25: computes `LSPA_table` using `LSPA_table=zeros(spin_mult,spin_mult)`.
- Lines 29: computes `LSPA_known` using `LSPA_known=T{n}*T{m}/norm(T{m},'fro')`.
- Lines 32: computes `RSPA_table` using `RSPA_table=zeros(spin_mult,spin_mult)`.
- Lines 36: computes `RSPA_known` using `RSPA_known=T{m}*T{n}/norm(T{m},'fro')`.
- Lines 39: computes `norm_diff_a` using `norm_diff_a=norm(LSPA_table-LSPA_known,'fro')`.
- Lines 40: computes `norm_diff_b` using `norm_diff_b=norm(RSPA_table-RSPA_known,'fro')`.
- Lines 41: computes `norm_scale` using `norm_scale=norm(T{m},'fro'); tol=sqrt(eps)`.

## Implementation structure

- Product action tests for irreducible spherical tensor
- operators and orthogonalised bosonic monomials.
- Random level counts
- Monomials and ISTs
- Multiplication tables
- IST test loop
- Left side product action
- Right side product action
- Accuracy tests
- BM test loop

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `randi()`, `irr_sph_ten()`, `boson_ortho()`, `ist_product_table()`, `bos_product_table()`, `ist_PTL()`, `ist_PTR()`, `bos_PTL()`, `bos_PTR()`.
