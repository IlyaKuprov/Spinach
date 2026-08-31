# examples/giant_spin/triple_dy_levels.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/triple_dy_levels.m`
- Signature: `triple_dy_levels()`
- Total lines: 124

## Purpose

Eight lowest energy levels as a function of the applied magnetic fi- eld in a triple Dy triangular complex -see Figure 12 in Ligand field parameters and g-tensor for the J=15/2 ground term were computed using the SINGLE_ANISO routine in MOLCAS. Calculation time: hours

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: Three J=15/2 dysprosium atoms; implemented by `sys.isotopes={'E16','E16','E16'}`.
- Lines 21-22: g-tensor eigenvalues; implemented by `g_eigs=[1.325781502 1.322640525 1.317917615]`.
- Lines 24-26: Spin-orbit corrections to the DD couplings; implemented by `sys.enable={'sodd'}`.
- Lines 28-29: g-tensor eigenvectors; implemented by `U3=[ -0.507708 0.520032 0.686877`.
- Lines 33-34: g-tensor matrix; implemented by `g3=U3'*diag(g_eigs)*U3`.
- Lines 36-37: Triangle arrangement; implemented by `R120=euler2dcm(0,0,-2*pi/3)`.
- Lines 40-41: Data supplied to Spinach; implemented by `inter.zeeman.matrix={g1,g2,g3}`.
- Lines 43-44: Coordinates; implemented by `inter.coordinates={[ 0.011931000 7.936473000 12.476785000]`.
- Lines 48-49: Exchange couplings (Spinach uses NMR convention); implemented by `J=icm2hz(0.0063)`.
- Lines 54-55: Rotate the ligand field into the molecular frame; implemented by `R=[-0.00244652355382 0.10885648156199 0.99405446578366`.
- Lines 60-61: Stevens coefficients; implemented by `Bkq{2}=[ 0.10326442519536E+01`.
- Lines 89-90: Convert to irreducible spherical tensors; implemented by `for k=2:2:6`.
- Lines 95-98: Supply to Spinach; implemented by `inter.giant.coeff={{[0 0 0],Bkq3{2},[0 0 0 0 0 0 0],Bkq3{4},[0 0 0 0 0 0 0 0 0 0 0],Bkq3{6}} {[0 0 0],Bkq3{2},[0 0 0 0 0 0 0],Bkq3{4},[0 0 0 0 0 0 0 0 0 0 0],Bkq3{6}} {[…`.
- Lines 103-104: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 107-108: This must be set to 1 Tesla; implemented by `sys.magnet=1.0`.
- Lines 110-111: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 114-115: Experiment parameters; implemented by `parameters.fields=[0 1]`.
- Lines 120-121: Run the simulation; implemented by `fieldscan_enlev(spin_system,parameters)`.

### Control flow inferred from the code

- Line 90: `for` loop over `k=2:2:6`.

### Key state/data transformations

- Lines 19: computes `sys.isotopes` using `sys.isotopes={'E16','E16','E16'}`.
- Lines 22: computes `g_eigs` using `g_eigs=[1.325781502 1.322640525 1.317917615]`.
- Lines 26: computes `sys.enable` using `sys.enable={'sodd'}`.
- Lines 29: computes `U3` using `U3=[ -0.507708 0.520032 0.686877`.
- Lines 34: computes `g3` using `g3=U3'*diag(g_eigs)*U3`.
- Lines 37: computes `R120` using `R120=euler2dcm(0,0,-2*pi/3)`.
- Lines 38: computes `g2` using `g2=R120*g3*R120'; g1=R120*g2*R120'`.
- Lines 41: computes `inter.zeeman.matrix` using `inter.zeeman.matrix={g1,g2,g3}`.
- Lines 44: computes `inter.coordinates` using `inter.coordinates={[ 0.011931000 7.936473000 12.476785000]`.
- Lines 49: computes `J` using `J=icm2hz(0.0063)`.
- Lines 50: computes `inter.coupling.scalar` using `inter.coupling.scalar={0 J J`.
- Lines 55: computes `R` using `R=[-0.00244652355382 0.10885648156199 0.99405446578366`.
- Lines 58: computes `[alp,bet,gam]` using `[alp,bet,gam]=dcm2euler(R')`.
- Lines 61: computes `Bkq{2}` using `Bkq{2}=[ 0.10326442519536E+01`.
- Lines 66: computes `Bkq{4}` using `Bkq{4}=[ 0.11295055284297E-01`.
- Lines 75: computes `Bkq{6}` using `Bkq{6}=[-0.21293172047501E-03`.
- Lines 91: computes `Bkq{k}` using `Bkq{k}=icm2hz(Bkq{k})`.
- Lines 92: computes `Bkq3{k}` using `Bkq3{k}=wigner(k,alp,bet,gam)*stev2sph(k,Bkq{k})`.

## Implementation structure

- Eight lowest energy levels as a function of the applied magnetic fi-
- eld in a triple Dy triangular complex -see Figure 12 in
- Ligand field parameters and g-tensor for the J=15/2 ground term were
- computed using the SINGLE_ANISO routine in MOLCAS.
- Calculation time: hours
- Three J=15/2 dysprosium atoms
- g-tensor eigenvalues
- Spin-orbit corrections
- to the DD couplings
- g-tensor eigenvectors
- g-tensor matrix
- Triangle arrangement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `euler2dcm()`, `icm2hz()`, `dcm2euler()`, `wigner()`, `stev2sph()`, `create()`, `basis()`, `fieldscan_enlev()`.
