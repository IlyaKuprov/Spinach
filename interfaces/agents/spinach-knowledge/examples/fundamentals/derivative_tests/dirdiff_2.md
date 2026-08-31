# examples/fundamentals/derivative_tests/dirdiff_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/derivative_tests/dirdiff_2.m`
- Signature: `dirdiff_2()`
- Total lines: 56

## Purpose

Test of matrix exponential differentiation of second order Magnus product quadrature (trapdiff.m) with the result com- pared to the central finite difference derivative. General coherent + non-symmetric dissipative case is tested.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Formalisms to test; implemented by `formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'}`.
- Lines 14-15: Loop over formalisms; implemented by `for n=1:numel(formalisms)`.
- Lines 17-18: Get the Spinach object; implemented by `spin_system=dirdiff_test_system(formalisms{n})`.
- Lines 20-21: Left and right drift generators, dissipative; implemented by `Hd={randn(50)+1i*randn(50), randn(50)+1i*randn(50)}`.
- Lines 23-24: Control operator; implemented by `Hc=randn(50)+1i*randn(50)`.
- Lines 26-27: A reasonable time step estimate; implemented by `dt=mean([1/norm(Hd{1},2), 1/norm(Hd{2},2)])`.
- Lines 29-30: Reasonable controls; implemented by `cL=0.1*randn()*norm(Hd{1},2)`.
- Lines 33-34: Get analytical derivatives; implemented by `[DL_anl,DR_anl]=trapdiff(spin_system,Hd,Hc,dt,cL,cR)`.
- Lines 36-37: Get numerical derivatives; implemented by `dc=sqrt(eps('double'))`.
- Lines 46-48: Test analytical against finite difference; implemented by `if (norm(DL_anl-DL_num,2)<10*sqrt(eps('double')))&& (norm(DR_anl-DR_num,2)<10*sqrt(eps('double')))`.

### Control flow inferred from the code

- Line 15: `for` loop over `n=1:numel(formalisms)`.
- Line 47: conditional branch on `(norm(DL_anl-DL_num,2)<10*sqrt(eps('double')))&&`.

### Key state/data transformations

- Lines 12: computes `formalisms` using `formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'}`.
- Lines 18: computes `spin_system` using `spin_system=dirdiff_test_system(formalisms{n})`.
- Lines 21: computes `Hd` using `Hd={randn(50)+1i*randn(50), randn(50)+1i*randn(50)}`.
- Lines 24: computes `Hc` using `Hc=randn(50)+1i*randn(50)`.
- Lines 27: computes `dt` using `dt=mean([1/norm(Hd{1},2), 1/norm(Hd{2},2)])`.
- Lines 30: computes `cL` using `cL=0.1*randn()*norm(Hd{1},2)`.
- Lines 31: computes `cR` using `cR=0.1*randn()*norm(Hd{2},2)`.
- Lines 34: computes `[DL_anl,DR_anl]` using `[DL_anl,DR_anl]=trapdiff(spin_system,Hd,Hc,dt,cL,cR)`.
- Lines 37: computes `dc` using `dc=sqrt(eps('double'))`.
- Lines 38: computes `H_dir_L` using `H_dir_L=(1/2)*Hc+1i*dt*(sqrt(3)/12)*(Hc*Hd{2}-Hd{2}*Hc)`.
- Lines 39: computes `H_dir_R` using `H_dir_R=(1/2)*Hc+1i*dt*(sqrt(3)/12)*(Hd{1}*Hc-Hc*Hd{1})`.
- Lines 40: computes `H` using `H=(Hd{1}+Hd{2})/2+cL*H_dir_L+cR*H_dir_R`.
- Lines 41-42: computes `DL_num` using `DL_num=(expm(-1i*dt*(H+dc*H_dir_L))- expm(-1i*dt*(H-dc*H_dir_L)))/(2*dc)`.
- Lines 43-44: computes `DR_num` using `DR_num=(expm(-1i*dt*(H+dc*H_dir_R))- expm(-1i*dt*(H-dc*H_dir_R)))/(2*dc)`.

## Implementation structure

- Test of matrix exponential differentiation of second order
- Magnus product quadrature (trapdiff.m) with the result com-
- pared to the central finite difference derivative. General
- coherent + non-symmetric dissipative case is tested.
- Formalisms to test
- Loop over formalisms
- Get the Spinach object
- Left and right drift generators, dissipative
- Control operator
- A reasonable time step estimate
- Reasonable controls
- Get analytical derivatives

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `dirdiff_test_system()`, `trapdiff()`, `eps()`.
