# examples/fundamentals/tensor_structures/polyadic_test_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/tensor_structures/polyadic_test_1.m`
- Signature: `polyadic_test_1()`
- Total lines: 120

## Purpose

Unit tests for the polyadic object.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 7-8: Get random complex matrices; implemented by `a=randn(7,7)+1i*randn(7,7)`.
- Lines 12-13: Get a random complex vector; implemented by `v=randn(7*9,1)+1i*randn(7*9,1)`.
- Lines 15-16: Normalise everything; implemented by `a=a/norm(a,2); b=b/norm(b,2)`.
- Lines 19-20: Create the polyadic; implemented by `P=polyadic({{polyadic({{a,b}}),c},{a b c}})`.
- Lines 22-23: Add prefixes and suffixes; implemented by `p1=randn(size(P))+1i*randn(size(P)); p1=p1/norm(p1,2)`.
- Lines 30-31: Reference matrix; implemented by `M=2*p1*p2*kron(a,kron(b,c))*s1*s2`.
- Lines 33-35: Create-inflate test; implemented by `if (norm(M-inflate(P),1)<1e-14)&& (norm(M-full(P),1)<1e-14)`.
- Lines 41-42: Matrix-vector test; implemented by `Pv=P*v; vP=v'*P; Mv=M*v; vM=v'*M`.
- Lines 50-52: Addition test 1; implemented by `if (norm(inflate(P-P),1)<1e-14)&& (norm(full(P-P),1)<1e-14)`.
- Lines 58-60: Addition test 2; implemented by `if (norm(inflate(2*P+P*3)-5*M,1)<1e-14)&& (norm(full(2*P+P*3)-5*M,1)<1e-14)`.
- Lines 66-67: Conjugate-transpose test 1; implemented by `if norm(P*v-(v'*P')')<1e-14`.
- Lines 73-74: Conjugate-transpose test 2; implemented by `if norm((inflate(P)')*v-inflate(P')*v)<1e-14`.
- Lines 80-81: Size test; implemented by `if norm(size(inflate(P))-size(P))<1e-14`.
- Lines 87-88: Kron test 1; implemented by `K=randn(2,2)+1i*randn(2,2); K=K/norm(K,2)`.
- Lines 96-100: Kron test 2; implemented by `K=polyadic({{randn(2,2)+1i*randn(2,2), randn(2,2)+1i*randn(2,2)}, {randn(2,2)+1i*randn(2,2), randn(2,2)+1i*randn(2,2)}}); K=(1/norm(full(K),2))*K`.
- Lines 108-109: Step test; implemented by `spin_system=bootstrap('hush')`.

### Control flow inferred from the code

- Line 34: conditional branch on `(norm(M-inflate(P),1)<1e-14)&&`.
- Line 43: conditional branch on `(norm(Pv-Mv,1)<1e-14)&&`.
- Line 51: conditional branch on `(norm(inflate(P-P),1)<1e-14)&&`.
- Line 59: conditional branch on `(norm(inflate(2*P+P*3)-5*M,1)<1e-14)&&`.
- Line 67: conditional branch on `norm(P*v-(v'*P')')<1e-14`.
- Line 74: conditional branch on `norm((inflate(P)')*v-inflate(P')*v)<1e-14`.
- Line 81: conditional branch on `norm(size(inflate(P))-size(P))<1e-14`.
- Line 89: conditional branch on `(norm(kron(inflate(P),K)-inflate(kron(P,K)),1)<1e-14)&&`.
- Line 101: conditional branch on `(norm(kron(full(K),inflate(P))-inflate(kron(K,P)),1)<1e-14)&&`.
- Line 113: conditional branch on `norm(v1-v2)<1e-12`.

### Key state/data transformations

- Lines 8: computes `a` using `a=randn(7,7)+1i*randn(7,7)`.
- Lines 9: computes `b` using `b=randn(1,1)+1i*randn(1,1)`.
- Lines 10: computes `c` using `c=sprandn(9,9,0.1)+1i*sprandn(9,9,0.1)`.
- Lines 13: computes `v` using `v=randn(7*9,1)+1i*randn(7*9,1)`.
- Lines 20: computes `P` using `P=polyadic({{polyadic({{a,b}}),c},{a b c}})`.
- Lines 23: computes `p1` using `p1=randn(size(P))+1i*randn(size(P)); p1=p1/norm(p1,2)`.
- Lines 24: computes `p2` using `p2=randn(size(P))+1i*randn(size(P)); p2=p2/norm(p2,2)`.
- Lines 25: computes `s1` using `s1=randn(size(P))+1i*randn(size(P)); s1=s1/norm(s1,2)`.
- Lines 26: computes `s2` using `s2=randn(size(P))+1i*randn(size(P)); s2=s2/norm(s2,2)`.
- Lines 31: computes `M` using `M=2*p1*p2*kron(a,kron(b,c))*s1*s2`.
- Lines 42: computes `Pv` using `Pv=P*v; vP=v'*P; Mv=M*v; vM=v'*M`.
- Lines 88: computes `K` using `K=randn(2,2)+1i*randn(2,2); K=K/norm(K,2)`.
- Lines 99-100: computes `{randn(2,2)+1i*randn(2,2), randn(2,2)+1i*randn(2,2)}}); K` using `{randn(2,2)+1i*randn(2,2), randn(2,2)+1i*randn(2,2)}}); K=(1/norm(full(K),2))*K`.
- Lines 100: computes `randn(2,2)+1i*randn(2,2)}}); K` using `randn(2,2)+1i*randn(2,2)}}); K=(1/norm(full(K),2))*K`.
- Lines 109: computes `spin_system` using `spin_system=bootstrap('hush')`.
- Lines 110: computes `dt` using `dt=1/norm(M,1)`.
- Lines 111: computes `v1` using `v1=step(spin_system,P,v,dt)`.
- Lines 112: computes `v2` using `v2=step(spin_system,M,v,dt)`.

## Implementation structure

- Unit tests for the polyadic object.
- Get random complex matrices
- Get a random complex vector
- Normalise everything
- Create the polyadic
- Add prefixes and suffixes
- Reference matrix
- Create-inflate test
- Matrix-vector test
- Addition test 1
- Addition test 2
- Conjugate-transpose test 1

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sprandn()`, `polyadic()`, `prefix()`, `suffix()`, `inflate()`, `bootstrap()`, `step()`.
