# kernel/conventions/transforms/dcm2wigner.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/dcm2wigner.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=dcm2wigner.m)

## Conversion

For a real 3x3 directional cosine matrix dcm, the function constructs the rank-2 Wigner matrix D. Its rows and columns are ordered by magnetic index m = 2, 1, 0, -1, -2. The documented action is v_out = D*v_in, where v_in contains irreducible spherical tensor coefficients in the same order: T(2,2), T(2,1), T(2,0), T(2,-1), T(2,-2). The source assigns no physical units to D and does not specify whether dcm represents an active or passive rotation; no such convention is inferred here.

The implementation defines complex scalars A and B from the matrix entries, then Z:

~~~text
A = sqrt(0.5*(dcm(1,1) + i*dcm(1,2) - i*dcm(2,1) + dcm(2,2)))
B = sqrt(0.5*(-dcm(1,1) + i*dcm(1,2) + i*dcm(2,1) + dcm(2,2)))
Z = A*conj(A) - B*conj(B)
~~~

Here conj denotes scalar complex conjugation (MATLAB apostrophe in the source). With these coefficients, the returned matrix is:

~~~text
D = [ A^4                    2*A^3*B                 sqrt(6)*A^2*B^2      2*A*B^3                    B^4
     -2*A^3*conj(B)          A^2*(2*Z-1)             sqrt(6)*A*B*Z        B^2*(2*Z+1)                2*conj(A)*B^3
      sqrt(6)*A^2*conj(B)^2 -sqrt(6)*A*conj(B)*Z    0.5*(3*Z^2-1)        sqrt(6)*conj(A)*B*Z       sqrt(6)*conj(A)^2*B^2
     -2*A*conj(B)^3          conj(B)^2*(2*Z+1)      -sqrt(6)*conj(A)*conj(B)*Z  conj(A)^2*(2*Z-1)  2*conj(A)^3*B
      conj(B)^4              -2*conj(A)*conj(B)^3  sqrt(6)*conj(A)^2*conj(B)^2 -2*conj(A)^3*conj(B)  conj(A)^4 ];
~~~

## Input and checks

The input must be a real numeric 3x3 matrix. The validator warns if norm(dcm'*dcm-eye(3),1)>1e-6 or abs(det(dcm)-1)>1e-6; either quantity above 1e-2 raises an error. The conversion also checks amplitude and phase consistency at 1e-6; a phase mismatch first changes the sign of A, and a remaining mismatch raises an error. The source has no separate finite-value check.
