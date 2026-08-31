# examples/fundamentals/operator_tests/commutation_6.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/operator_tests/commutation_6.m`
- Signature: `commutation_6()`
- Total lines: 91

## Purpose

Commutation and product relations for Pauli and central transition operators.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Accuracy threshold; implemented by `tol=1e-10`.
- Lines 11-12: Test SU(2) relations for Pauli operators; implemented by `for mult=[2 3 4 7]`.
- Lines 14-15: Build Pauli operators; implemented by `S=pauli(mult)`.
- Lines 17-18: Check angular momentum commutators; implemented by `err_xy=norm(comm(S.x,S.y)-1i*S.z,'fro')`.
- Lines 22-23: Check ladder operator definitions; implemented by `err_p=norm(S.p-(S.x+1i*S.y),'fro')`.
- Lines 26-27: Report Pauli operator failures; implemented by `if max([err_xy err_yz err_zx err_p err_m])>tol`.
- Lines 34-35: Test SU(2) relations for central transition operators; implemented by `for mult=[4 6 8]`.
- Lines 37-38: Build central transition operators; implemented by `CTx=centrans(mult,'x')`.
- Lines 44-45: Build central transition support projector; implemented by `P=zeros(mult,mult)`.
- Lines 49-50: Check commutators and product identities; implemented by `err_comm=norm(comm(CTx,CTy)-1i*CTz,'fro')`.
- Lines 56-57: Report central transition failures; implemented by `if max([err_comm err_lad_1 err_lad_2 err_lad_3 err_prod])>tol`.
- Lines 64-65: Test central transition IST expansions; implemented by `for mult=[4 6 8]`.
- Lines 68-69: Build central transition operator; implemented by `C=centrans(mult,type{1})`.
- Lines 71-72: Obtain IST expansion.; implemented by `[states,coeffs]=ct2ist(mult,type{1})`.
- Lines 75-76: Reconstruct the operator from IST basis terms; implemented by `C_rec=zeros(mult,mult,'like',1i)`.
- Lines 81-82: Report IST expansion failures; implemented by `if norm(C-C_rec,'fro')>tol`.

### Control flow inferred from the code

- Line 12: `for` loop over `mult=[2 3 4 7]`.
- Line 27: conditional branch on `max([err_xy err_yz err_zx err_p err_m])>tol`.
- Line 35: `for` loop over `mult=[4 6 8]`.
- Line 57: conditional branch on `max([err_comm err_lad_1 err_lad_2 err_lad_3 err_prod])>tol`.
- Line 65: `for` loop over `mult=[4 6 8]`.
- Line 66: `for` loop over `type={'x','y','z','+','-'}`.
- Line 77: `for` loop over `n=1:numel(states)`.
- Line 82: conditional branch on `norm(C-C_rec,'fro')>tol`.

### Key state/data transformations

- Lines 9: computes `tol` using `tol=1e-10`.
- Lines 15: computes `S` using `S=pauli(mult)`.
- Lines 18: computes `err_xy` using `err_xy=norm(comm(S.x,S.y)-1i*S.z,'fro')`.
- Lines 19: computes `err_yz` using `err_yz=norm(comm(S.y,S.z)-1i*S.x,'fro')`.
- Lines 20: computes `err_zx` using `err_zx=norm(comm(S.z,S.x)-1i*S.y,'fro')`.
- Lines 23: computes `err_p` using `err_p=norm(S.p-(S.x+1i*S.y),'fro')`.
- Lines 24: computes `err_m` using `err_m=norm(S.m-(S.x-1i*S.y),'fro')`.
- Lines 38: computes `CTx` using `CTx=centrans(mult,'x')`.
- Lines 39: computes `CTy` using `CTy=centrans(mult,'y')`.
- Lines 40: computes `CTz` using `CTz=centrans(mult,'z')`.
- Lines 41: computes `CTp` using `CTp=centrans(mult,'+')`.
- Lines 42: computes `CTm` using `CTm=centrans(mult,'-')`.
- Lines 45: computes `P` using `P=zeros(mult,mult)`.
- Lines 46: computes `P(mult/2,mult/2)` using `P(mult/2,mult/2)=1`.
- Lines 47: computes `P(mult/2+1,mult/2+1)` using `P(mult/2+1,mult/2+1)=1`.
- Lines 50: computes `err_comm` using `err_comm=norm(comm(CTx,CTy)-1i*CTz,'fro')`.
- Lines 51: computes `err_lad_1` using `err_lad_1=norm(comm(CTz,CTp)-CTp,'fro')`.
- Lines 52: computes `err_lad_2` using `err_lad_2=norm(comm(CTz,CTm)+CTm,'fro')`.

## Implementation structure

- Commutation and product relations for Pauli and
- central transition operators.
- Accuracy threshold
- Test SU(2) relations for Pauli operators
- Build Pauli operators
- Check angular momentum commutators
- Check ladder operator definitions
- Report Pauli operator failures
- Test SU(2) relations for central transition operators
- Build central transition operators
- Build central transition support projector
- Check commutators and product identities

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `pauli()`, `comm()`, `centrans()`, `ct2ist()`, `irr_sph_ten()`, `coeffs()`, `states()`.
