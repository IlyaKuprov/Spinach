# examples/fundamentals/operator_tests/commutation_7.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/operator_tests/commutation_7.m`
- Signature: `commutation_7()`
- Total lines: 126

## Purpose

Expansion relations for operator basis transforms.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 7-8: Accuracy threshold; implemented by `tol=1e-10`.
- Lines 10-11: Test IST expansion of Zeeman level projectors; implemented by `mult=5; T=irr_sph_ten(mult)`.
- Lines 14-15: Build known Zeeman projector; implemented by `P=zeros(mult,mult)`.
- Lines 18-19: Obtain IST expansion.; implemented by `[states,coeffs]=enlev2ist(mult,lvl_num,'S')`.
- Lines 21-22: Reconstruct operator from IST terms; implemented by `P_rec=zeros(mult,mult,'like',1i)`.
- Lines 27-28: Report IST expansion failures; implemented by `if norm(P-P_rec,'fro')>tol`.
- Lines 35-36: Test BM expansion of bosonic level projectors; implemented by `nlevels=6; B=boson_mono(nlevels)`.
- Lines 39-40: Build known bosonic projector; implemented by `P=zeros(nlevels,nlevels)`.
- Lines 43-44: Obtain BM expansion.; implemented by `[states,coeffs]=enlev2bm(nlevels,lvl_num)`.
- Lines 46-47: Reconstruct operator from BM terms; implemented by `P_rec=zeros(nlevels,nlevels,'like',1i)`.
- Lines 52-53: Report BM expansion failures; implemented by `if norm(P-P_rec,'fro')>tol`.
- Lines 60-61: Test bosonic product to IST conversion; implemented by `W=weyl(nlevels); T_bos=irr_sph_ten(nlevels)`.
- Lines 65-66: Build operator product explicitly; implemented by `P=speye(nlevels,nlevels)`.
- Lines 75-76: Obtain IST expansion; implemented by `[states,coeffs]=bos2ist(prod_list{n},nlevels)`.
- Lines 78-79: Reconstruct operator from IST terms; implemented by `P_rec=zeros(nlevels,nlevels,'like',1i)`.
- Lines 84-85: Report bosonic IST conversion failures; implemented by `if norm(P-P_rec,'fro')>tol`.
- Lines 92-93: Test matrix expansion into single-transition basis; implemented by `dim=7; A=randn(dim)+1i*randn(dim)`.
- Lines 97-98: Add transition basis term; implemented by `A_rec=A_rec+hdot(E{n},A)*E{n}`.

### Control flow inferred from the code

- Line 12: `for` loop over `lvl_num=1:mult`.
- Line 23: `for` loop over `n=1:numel(states)`.
- Line 28: conditional branch on `norm(P-P_rec,'fro')>tol`.
- Line 37: `for` loop over `lvl_num=1:nlevels`.
- Line 48: `for` loop over `n=1:numel(states)`.
- Line 53: conditional branch on `norm(P-P_rec,'fro')>tol`.
- Line 63: `for` loop over `n=1:numel(prod_list)`.
- Line 67: `for` loop over `k=1:numel(prod_list{n})`.
- Line 68: conditional branch on `prod_list{n}(k)=='C'`.
- Line 80: `for` loop over `k=1:numel(states)`.
- Line 85: conditional branch on `norm(P-P_rec,'fro')>tol`.
- Line 95: `for` loop over `n=1:numel(E)`.
- Line 103: conditional branch on `norm(A-A_rec,'fro')>tol`.
- Line 112: `for` loop over `n=1:numel(states)`.

### Key state/data transformations

- Lines 8: computes `tol` using `tol=1e-10`.
- Lines 11: computes `mult` using `mult=5; T=irr_sph_ten(mult)`.
- Lines 15: computes `P` using `P=zeros(mult,mult)`.
- Lines 16: computes `P(mult-lvl_num+1,mult-lvl_num+1)` using `P(mult-lvl_num+1,mult-lvl_num+1)=1`.
- Lines 19: computes `[states,coeffs]` using `[states,coeffs]=enlev2ist(mult,lvl_num,'S')`.
- Lines 22: computes `P_rec` using `P_rec=zeros(mult,mult,'like',1i)`.
- Lines 36: computes `nlevels` using `nlevels=6; B=boson_mono(nlevels)`.
- Lines 41: computes `P(lvl_num,lvl_num)` using `P(lvl_num,lvl_num)=1`.
- Lines 61: computes `W` using `W=weyl(nlevels); T_bos=irr_sph_ten(nlevels)`.
- Lines 62: computes `prod_list` using `prod_list={'CA','ACCA','CCAAA'}`.
- Lines 93: computes `dim` using `dim=7; A=randn(dim)+1i*randn(dim)`.
- Lines 94: computes `E` using `E=sin_tran(dim); A_rec=zeros(dim,dim,'like',1i)`.
- Lines 98: computes `A_rec` using `A_rec=A_rec+hdot(E{n},A)*E{n}`.
- Lines 109: computes `A` using `A=randn(nlevels)+1i*randn(nlevels)`.

## Implementation structure

- Expansion relations for operator basis transforms.
- Accuracy threshold
- Test IST expansion of Zeeman level projectors
- Build known Zeeman projector
- Obtain IST expansion.
- Reconstruct operator from IST terms
- Report IST expansion failures
- Test BM expansion of bosonic level projectors
- Build known bosonic projector
- Obtain BM expansion.
- Reconstruct operator from BM terms
- Report BM expansion failures

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `irr_sph_ten()`, `enlev2ist()`, `coeffs()`, `states()`, `boson_mono()`, `enlev2bm()`, `weyl()`, `speye()`, `bos2ist()`, `sin_tran()`, `hdot()`, `oper2bm()`.
