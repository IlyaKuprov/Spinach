# examples/giant_spin/triple_dy_magn.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/triple_dy_magn.m`
- Signature: `triple_dy_magn()`
- Total lines: 133

## Purpose

Simulation of a finite-speed magnetic field sweep experiment for a single crystal of a triple-Dy triangular complex in a micro-SQUID, see Figure S24 in the Supplementary Information of Ligand field parameters and g-tensor for the J=15/2 ground term were computed using the SINGLE_ANISO routine in MOLCAS. Calculation time: hours

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Three J=15/2 dysprosium atoms; implemented by `sys.isotopes={'E16','E16','E16'}`.
- Lines 20-21: g-tensor eigenvalues; implemented by `g_eigs=[1.325781502 1.322640525 1.317917615]`.
- Lines 23-25: Spin-orbit corrections to the DD couplings; implemented by `sys.enable={'sodd'}`.
- Lines 27-28: g-tensor eigenvectors; implemented by `U3=[ -0.507708 0.520032 0.686877`.
- Lines 32-33: g-tensor matrix; implemented by `g3=U3'*diag(g_eigs)*U3`.
- Lines 35-36: Triangle arrangement; implemented by `R120=euler2dcm(0,0,-2*pi/3)`.
- Lines 39-40: Data supplied to Spinach; implemented by `inter.zeeman.matrix={g1,g2,g3}`.
- Lines 42-43: Coordinates; implemented by `inter.coordinates={[ 0.011931000 7.936473000 12.476785000]`.
- Lines 47-48: Exchange couplings (Spinach uses NMR convention); implemented by `J=icm2hz(0.0063)`.
- Lines 53-54: Rotate the ligand field into the molecular frame; implemented by `R=[-0.00244652355382 0.10885648156199 0.99405446578366`.
- Lines 59-60: Stevens coefficients; implemented by `Bkq{2}=[ 0.10326442519536E+01`.
- Lines 88-89: Convert to irreducible spherical tensors; implemented by `for k=2:2:6`.
- Lines 94-97: Supply to Spinach; implemented by `inter.giant.coeff={{[0 0 0],Bkq3{2},[0 0 0 0 0 0 0],Bkq3{4},[0 0 0 0 0 0 0 0 0 0 0],Bkq3{6}} {[0 0 0],Bkq3{2},[0 0 0 0 0 0 0],Bkq3{4},[0 0 0 0 0 0 0 0 0 0 0],Bkq3{6}} {[…`.
- Lines 102-103: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 106-107: Temperature in Kelvin; implemented by `inter.temperature=0.03`.
- Lines 109-110: This must be set to 1 Tesla; implemented by `sys.magnet=1.0`.
- Lines 112-113: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 116-117: Experiment parameters; implemented by `parameters.fields=[0 1]`.

### Control flow inferred from the code

- Line 89: `for` loop over `k=2:2:6`.

### Key state/data transformations

- Lines 18: computes `sys.isotopes` using `sys.isotopes={'E16','E16','E16'}`.
- Lines 21: computes `g_eigs` using `g_eigs=[1.325781502 1.322640525 1.317917615]`.
- Lines 25: computes `sys.enable` using `sys.enable={'sodd'}`.
- Lines 28: computes `U3` using `U3=[ -0.507708 0.520032 0.686877`.
- Lines 33: computes `g3` using `g3=U3'*diag(g_eigs)*U3`.
- Lines 36: computes `R120` using `R120=euler2dcm(0,0,-2*pi/3)`.
- Lines 37: computes `g2` using `g2=R120*g3*R120'; g1=R120*g2*R120'`.
- Lines 40: computes `inter.zeeman.matrix` using `inter.zeeman.matrix={g1,g2,g3}`.
- Lines 43: computes `inter.coordinates` using `inter.coordinates={[ 0.011931000 7.936473000 12.476785000]`.
- Lines 48: computes `J` using `J=icm2hz(0.0063)`.
- Lines 49: computes `inter.coupling.scalar` using `inter.coupling.scalar={0 J J`.
- Lines 54: computes `R` using `R=[-0.00244652355382 0.10885648156199 0.99405446578366`.
- Lines 57: computes `[alp,bet,gam]` using `[alp,bet,gam]=dcm2euler(R')`.
- Lines 60: computes `Bkq{2}` using `Bkq{2}=[ 0.10326442519536E+01`.
- Lines 65: computes `Bkq{4}` using `Bkq{4}=[ 0.11295055284297E-01`.
- Lines 74: computes `Bkq{6}` using `Bkq{6}=[-0.21293172047501E-03`.
- Lines 90: computes `Bkq{k}` using `Bkq{k}=icm2hz(Bkq{k})`.
- Lines 91: computes `Bkq3{k}` using `Bkq3{k}=wigner(k,alp,bet,gam)*stev2sph(k,Bkq{k})`.

## Implementation structure

- Simulation of a finite-speed magnetic field sweep experiment for a
- single crystal of a triple-Dy triangular complex in a micro-SQUID,
- see Figure S24 in the Supplementary Information of
- Ligand field parameters and g-tensor for the J=15/2 ground term were
- computed using the SINGLE_ANISO routine in MOLCAS.
- Calculation time: hours
- Three J=15/2 dysprosium atoms
- g-tensor eigenvalues
- Spin-orbit corrections
- to the DD couplings
- g-tensor eigenvectors
- g-tensor matrix

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `euler2dcm()`, `icm2hz()`, `dcm2euler()`, `wigner()`, `stev2sph()`, `create()`, `basis()`, `fieldscan_magn()`, `kfigure()`, `kxlabel()`, `kylabel()`.
