# examples/fundamentals/state_tests/state_consistency_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_tests/state_consistency_3.m`
- Signature: `state_consistency_3()`
- Total lines: 81

## Purpose

Deuterium pair singlet, triplet, and quintet state internal consistency test.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: A pair of deuteria; implemented by `sys.magnet=0`.
- Lines 13-14: Hilbert space; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 17-18: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 21-22: Ortho-deuterium states from Spinach; implemented by `[S,T,Q]=deut_pair(spin_system,1,2)`.
- Lines 24-26: Component vectors in Hilbert space as per Eq 1; implemented by `alp=[1; 0; 0]; bet=[0; 1; 0]; gam=[0; 0; 1]`.
- Lines 37-38: Test singlet state; implemented by `if norm(S-S0*S0',1)>1e-6`.
- Lines 42-45: Test triplet states; implemented by `if (norm(T{1}-Tp*Tp',1)>1e-6)|| (norm(T{2}-T0*T0',1)>1e-6)|| (norm(T{3}-Tm*Tm',1)>1e-6)`.
- Lines 49-54: Test quintet states; implemented by `if (norm(Q{1}-Qpp*Qpp',1)>1e-6)|| (norm(Q{2}-Qp*Qp',1)>1e-6)|| (norm(Q{3}-Q0*Q0',1)>1e-6)|| (norm(Q{4}-Qm*Qm',1)>1e-6)|| (norm(Q{5}-Qmm*Qmm',1)>1e-6)`.
- Lines 58-59: Move to Liouville space; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 69-70: Unit state from Spinach; implemented by `U=state(spin_system,{'E','E'},{1 2})`.
- Lines 72-73: Summation into the unit state test; implemented by `if norm(S+T{1}+T{2}+T{3}+Q{1}+Q{2}+Q{3}+Q{4}+Q{5}-U,1)>1e-6`.
- Lines 77-78: Report success; implemented by `disp('State construction test PASSED.')`.

### Control flow inferred from the code

- Line 38: conditional branch on `norm(S-S0*S0',1)>1e-6`.
- Line 43: conditional branch on `(norm(T{1}-Tp*Tp',1)>1e-6)||`.
- Line 50: conditional branch on `(norm(Q{1}-Qpp*Qpp',1)>1e-6)||`.
- Line 73: conditional branch on `norm(S+T{1}+T{2}+T{3}+Q{1}+Q{2}+Q{3}+Q{4}+Q{5}-U,1)>1e-6`.

### Key state/data transformations

- Lines 9: computes `sys.magnet` using `sys.magnet=0`.
- Lines 10: computes `sys.isotopes` using `sys.isotopes={'2H','2H'}`.
- Lines 11: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0 0}`.
- Lines 14: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 15: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 18: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 22: computes `[S,T,Q]` using `[S,T,Q]=deut_pair(spin_system,1,2)`.
- Lines 26: computes `alp` using `alp=[1; 0; 0]; bet=[0; 1; 0]; gam=[0; 0; 1]`.
- Lines 27: computes `S0` using `S0=(1/sqrt(3))*(kron(alp,gam)-kron(bet,bet)+kron(gam,alp))`.
- Lines 28: computes `Tp` using `Tp=(1/sqrt(2))*(kron(alp,bet)-kron(bet,alp))`.
- Lines 29: computes `T0` using `T0=(1/sqrt(2))*(kron(alp,gam)-kron(gam,alp))`.
- Lines 30: computes `Tm` using `Tm=(1/sqrt(2))*(kron(bet,gam)-kron(gam,bet))`.
- Lines 31: computes `Qpp` using `Qpp=kron(alp,alp)`.
- Lines 32: computes `Qp` using `Qp=(1/sqrt(2))*(kron(alp,bet)+kron(bet,alp))`.
- Lines 33: computes `Q0` using `Q0=(1/sqrt(6))*(kron(alp,gam)+2*kron(bet,bet)+kron(gam,alp))`.
- Lines 34: computes `Qm` using `Qm=(1/sqrt(2))*(kron(bet,gam)+kron(gam,bet))`.
- Lines 35: computes `Qmm` using `Qmm=kron(gam,gam)`.
- Lines 70: computes `U` using `U=state(spin_system,{'E','E'},{1 2})`.

## Implementation structure

- Deuterium pair singlet, triplet, and quintet state
- internal consistency test.
- A pair of deuteria
- Hilbert space
- Spinach housekeeping
- Ortho-deuterium states from Spinach
- Component vectors in Hilbert space as per Eq 1
- in https://doi.org/10.1016/S0009-2614(98)00784-2
- Test singlet state
- Test triplet states
- Test quintet states
- Move to Liouville space

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `deut_pair()`, `state()`.
