# examples/giant_spin/nuclear_relaxation_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/nuclear_relaxation_1.m`
- Signature: `nuclear_relaxation_1()`
- Total lines: 148

## Purpose

Nuclear relaxation rates using the adiabatic elimination method for a rapidly relaxing Dy(III) ion with a user-specified ZFS. Calculation time: minutes

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Magnetic field; implemented by `sys.magnet=14.1`.
- Lines 14-15: Dy(III) ion and a proton; implemented by `sys.isotopes={'E16','1H'}`.
- Lines 17-18: Electron g-tensor; implemented by `D=[ 1.325781 1.322640 1.317917]`.
- Lines 24-26: Spin-orbit corrections to the DD couplings; implemented by `sys.enable={'sodd'}`.
- Lines 28-29: Nuclear shift tensor; implemented by `inter.zeeman.matrix{2}=zeros(3,3)`.
- Lines 31-32: Rotate the ligand field into the molecular frame; implemented by `R=[-0.00244652355382 0.10885648156199 0.99405446578366`.
- Lines 37-38: Liza -this needs more decimal places; implemented by `rkd=[ -0.2821 0.7292 0.6234`.
- Lines 43-44: Ligand field parameters (MOLCAS); implemented by `Bkq{2}=[ 0.10326442519536E+01`.
- Lines 72-73: Convert to irreducible spherical tensors; implemented by `for k=2:2:6`.
- Lines 79-82: Supply to Spinach; implemented by `inter.giant.coeff={{[0 0 0],SBkq{2}, [0 0 0 0 0 0 0],SBkq{4}, [0 0 0 0 0 0 0 0 0 0 0],SBkq{6}},{}}`.
- Lines 86-87: Cartesian coordinates; implemented by `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 90-91: Separate spin relaxation times; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 97-98: Formalism specification; implemented by `bas.formalism='sphten-liouv'`.
- Lines 101-102: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 105-106: Slow and fast subspace indices; implemented by `slow_idx=1:4`.
- Lines 109-110: Hamiltonian; implemented by `[I,Q]=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 112-113: Non-interacting relaxation; implemented by `R_ni=relaxation(spin_system)`.
- Lines 115-118: Get spherical averaging grid as a structure; implemented by `sph_grid=load([spin_system.sys.root_dir '/kernel/grids/' 'leb_2ang_rank_11.mat'],'alphas','betas', 'gammas','weights')`.

### Control flow inferred from the code

- Line 73: `for` loop over `k=2:2:6`.
- Line 124: `parfor` loop over `n=1:numel(weights)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'E16','1H'}`.
- Lines 18: computes `D` using `D=[ 1.325781 1.322640 1.317917]`.
- Lines 19: computes `V` using `V=[-0.507708 0.520032 0.686877`.
- Lines 22: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=V'*diag(D)*V`.
- Lines 26: computes `sys.enable` using `sys.enable={'sodd'}`.
- Lines 29: computes `inter.zeeman.matrix{2}` using `inter.zeeman.matrix{2}=zeros(3,3)`.
- Lines 32: computes `R` using `R=[-0.00244652355382 0.10885648156199 0.99405446578366`.
- Lines 35: computes `[alp,bet,gam]` using `[alp,bet,gam]=dcm2euler(R')`.
- Lines 38: computes `rkd` using `rkd=[ -0.2821 0.7292 0.6234`.
- Lines 41: computes `[a,b,g]` using `[a,b,g]=dcm2euler(rkd')`.
- Lines 44: computes `Bkq{2}` using `Bkq{2}=[ 0.10326442519536E+01`.
- Lines 49: computes `Bkq{4}` using `Bkq{4}=[ 0.11295055284297E-01`.
- Lines 58: computes `Bkq{6}` using `Bkq{6}=[-0.21293172047501E-03`.
- Lines 74: computes `Bkq{k}` using `Bkq{k}=icm2hz(Bkq{k})`.
- Lines 75-76: computes `SBkq{k}` using `SBkq{k}=wigner(k,a,b,g)* wigner(k,alp,bet,gam)*stev2sph(k,Bkq{k})`.
- Lines 80-82: computes `inter.giant.coeff` using `inter.giant.coeff={{[0 0 0],SBkq{2}, [0 0 0 0 0 0 0],SBkq{4}, [0 0 0 0 0 0 0 0 0 0 0],SBkq{6}},{}}`.
- Lines 83-84: computes `inter.giant.euler` using `inter.giant.euler={{[0 0 0],[0 0 0],[0 0 0], [0 0 0],[0 0 0],[0 0 0]},{}}`.

## Implementation structure

- Nuclear relaxation rates using the adiabatic elimination method
- for a rapidly relaxing Dy(III) ion with a user-specified ZFS.
- Calculation time: minutes
- Magnetic field
- Dy(III) ion and a proton
- Electron g-tensor
- Spin-orbit corrections
- to the DD couplings
- Nuclear shift tensor
- Rotate the ligand field into the molecular frame
- Liza -this needs more decimal places
- Ligand field parameters (MOLCAS)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `dcm2euler()`, `icm2hz()`, `wigner()`, `stev2sph()`, `create()`, `basis()`, `hamiltonian()`, `assume()`, `relaxation()`, `load()`, `orientation()`, `alphas()`, `betas()`, `gammas()`, `adelim()`, `weights()`.
