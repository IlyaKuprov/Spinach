# examples/giant_spin/triple_tb_eqmag_field.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/triple_tb_eqmag_field.m`
- Signature: `triple_tb_eqmag_field()`
- Total lines: 222

## Purpose

Simulation of the field dependence of the magnetisation of a triple- Tb triangular complex -see Figure S27 and S28 in the Supplementary Information of the following paper Ligand field parameters and g-tensor for the J=6 ground term were computed using the SINGLE_ANISO routine in MOLCAS. Calculation time: hours

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Three J=6 terbium atoms; implemented by `sys.isotopes={'E13','E13','E13'}`.
- Lines 20-21: g-tensor eigenvalues; implemented by `g1v=[1.497075749 1.495252923 1.481370349]`.
- Lines 25-27: Spin-orbit corrections to the DD couplings; implemented by `sys.enable={'sodd'}`.
- Lines 29-30: g-tensor eigenvalues (rows); implemented by `U1=[-0.847221 -0.510119 -0.148306`.
- Lines 40-41: g-tensor matrices; implemented by `g1=U1*diag(g1v)*U1'`.
- Lines 46-47: Tb atom coordinates; implemented by `inter.coordinates={[ 2.074806000 1.089908000 -0.432651000]`.
- Lines 51-52: Exchange couplings (Spinach uses NMR convention); implemented by `J=icm2hz(0.003)`.
- Lines 57-58: Ligand field rotations; implemented by `R1=[-0.84722146280957 -0.51011886371357 -0.14830555565586`.
- Lines 71-72: Stevens coefficients, Tb 1; implemented by `Bkq1{2}=[ 0.14292463619745E+00`.
- Lines 100-101: Stevens coefficients, Tb 2; implemented by `Bkq2{2}=[-0.51308162627882E+00`.
- Lines 129-130: Stevens coefficients, Tb 3; implemented by `Bkq3{2}=[ 0.60574413436613E+00`.
- Lines 158-159: Convert to irreducible spherical tensors; implemented by `for k=2:2:6`.
- Lines 168-171: Supply to Spinach; implemented by `inter.giant.coeff={{[0 0 0],Bkq1{2},[0 0 0 0 0 0 0],Bkq1{4},[0 0 0 0 0 0 0 0 0 0 0],Bkq1{6}} {[0 0 0],Bkq2{2},[0 0 0 0 0 0 0],Bkq2{4},[0 0 0 0 0 0 0 0 0 0 0],Bkq2{6}} {[…`.
- Lines 176-177: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 180-181: Spherical powder grid; implemented by `parameters.grid='leb_2ang_rank_11'`.
- Lines 183-184: Temperature in Kelvin; implemented by `inter.temperature=2.0`.
- Lines 186-188: Magnetic field range; implemented by `B0=[0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 2 3 4 5 6]`.
- Lines 190-191: Preallocate the answer; implemented by `Mz=nan(size(B0)); kfigure()`.

### Control flow inferred from the code

- Line 159: `for` loop over `k=2:2:6`.
- Line 194: `for` loop over `n=1:numel(B0)`.

### Key state/data transformations

- Lines 18: computes `sys.isotopes` using `sys.isotopes={'E13','E13','E13'}`.
- Lines 21: computes `g1v` using `g1v=[1.497075749 1.495252923 1.481370349]`.
- Lines 22: computes `g2v` using `g2v=[1.496540374 1.494559188 1.482686210]`.
- Lines 23: computes `g3v` using `g3v=[1.497265940 1.494866241 1.481858735]`.
- Lines 27: computes `sys.enable` using `sys.enable={'sodd'}`.
- Lines 30: computes `U1` using `U1=[-0.847221 -0.510119 -0.148306`.
- Lines 33: computes `U2` using `U2=[ 0.632422 -0.325205 0.703053`.
- Lines 36: computes `U3` using `U3=[ 0.038754 -0.801228 0.597102`.
- Lines 41: computes `g1` using `g1=U1*diag(g1v)*U1'`.
- Lines 42: computes `g2` using `g2=U2*diag(g2v)*U2'`.
- Lines 43: computes `g3` using `g3=U3*diag(g3v)*U3'`.
- Lines 44: computes `inter.zeeman.matrix` using `inter.zeeman.matrix={g1,g2,g3}`.
- Lines 47: computes `inter.coordinates` using `inter.coordinates={[ 2.074806000 1.089908000 -0.432651000]`.
- Lines 52: computes `J` using `J=icm2hz(0.003)`.
- Lines 53: computes `inter.coupling.scalar` using `inter.coupling.scalar={0 J J`.
- Lines 58: computes `R1` using `R1=[-0.84722146280957 -0.51011886371357 -0.14830555565586`.
- Lines 61: computes `R2` using `R2=[ 0.63242234003580 -0.32520539385769 0.70305293942172`.
- Lines 64: computes `R3` using `R3=[ 0.03875378043232 -0.80122835625537 0.59710239124837`.

## Implementation structure

- Simulation of the field dependence of the magnetisation of a triple-
- Tb triangular complex -see Figure S27 and S28 in the Supplementary
- Information of the following paper
- Ligand field parameters and g-tensor for the J=6 ground term were
- computed using the SINGLE_ANISO routine in MOLCAS.
- Calculation time: hours
- Three J=6 terbium atoms
- g-tensor eigenvalues
- Spin-orbit corrections
- to the DD couplings
- g-tensor eigenvalues (rows)
- g-tensor matrices

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `icm2hz()`, `dcm2euler()`, `wigner()`, `stev2sph()`, `nan()`, `kfigure()`, `create()`, `basis()`, `eqmag()`, `mag()`, `kxlabel()`, `kylabel()`, `load()`.
