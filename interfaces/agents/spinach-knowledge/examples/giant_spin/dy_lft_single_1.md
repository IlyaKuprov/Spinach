# examples/giant_spin/dy_lft_single_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/dy_lft_single_1.m`
- Signature: `dy_lft_single_1()`
- Total lines: 97

## Purpose

Reproduction of MOLCAS results with the Ligand Field Theory model for a single Dy(III) ion. Calculation time: seconds

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Magnetic field; implemented by `sys.magnet=0`.
- Lines 14-15: Single Dy ion; implemented by `sys.isotopes={'E16'}`.
- Lines 17-18: Real g-tensor; implemented by `D=[ 1.325781 1.322640 1.317917]`.
- Lines 24-25: Rotate the ligand field into the molecular frame; implemented by `R=[-0.00244652355382 0.10885648156199 0.99405446578366`.
- Lines 30-31: Liza -this needs more decimal places; implemented by `rkd=[ -0.2821 0.7292 0.6234`.
- Lines 36-37: Ligand field parameters (MOLCAS); implemented by `Bkq{2}=[ 0.10326442519536E+01`.
- Lines 65-66: Convert to irreducible spherical tensors; implemented by `for k=2:2:6`.
- Lines 72-75: Supply to Spinach; implemented by `inter.giant.coeff={{[0 0 0],SBkq{2}, [0 0 0 0 0 0 0],SBkq{4}, [0 0 0 0 0 0 0 0 0 0 0],SBkq{6}}}`.
- Lines 79-80: Formalism specification; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 83-84: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 87-88: Display Spinach answer; implemented by `[~,D]=eig(geffect(spin_system,[1 2]))`.
- Lines 92-93: Display MOLCAS answer; implemented by `disp(' '); disp('MOLCAS results:')`.

### Control flow inferred from the code

- Line 66: `for` loop over `k=2:2:6`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=0`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'E16'}`.
- Lines 18: computes `D` using `D=[ 1.325781 1.322640 1.317917]`.
- Lines 19: computes `V` using `V=[-0.507708 0.520032 0.686877`.
- Lines 22: computes `inter.zeeman.matrix` using `inter.zeeman.matrix={V'*diag(D)*V}`.
- Lines 25: computes `R` using `R=[-0.00244652355382 0.10885648156199 0.99405446578366`.
- Lines 28: computes `[alp,bet,gam]` using `[alp,bet,gam]=dcm2euler(R')`.
- Lines 31: computes `rkd` using `rkd=[ -0.2821 0.7292 0.6234`.
- Lines 34: computes `[a,b,g]` using `[a,b,g]=dcm2euler(rkd')`.
- Lines 37: computes `Bkq{2}` using `Bkq{2}=[ 0.10326442519536E+01`.
- Lines 42: computes `Bkq{4}` using `Bkq{4}=[ 0.11295055284297E-01`.
- Lines 51: computes `Bkq{6}` using `Bkq{6}=[-0.21293172047501E-03`.
- Lines 67: computes `Bkq{k}` using `Bkq{k}=icm2hz(Bkq{k})`.
- Lines 68-69: computes `SBkq{k}` using `SBkq{k}=wigner(k,a,b,g)* wigner(k,alp,bet,gam)*stev2sph(k,Bkq{k})`.
- Lines 73-75: computes `inter.giant.coeff` using `inter.giant.coeff={{[0 0 0],SBkq{2}, [0 0 0 0 0 0 0],SBkq{4}, [0 0 0 0 0 0 0 0 0 0 0],SBkq{6}}}`.
- Lines 76-77: computes `inter.giant.euler` using `inter.giant.euler={{[0 0 0],[0 0 0],[0 0 0], [0 0 0],[0 0 0],[0 0 0]}}`.
- Lines 80: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 81: computes `bas.approximation` using `bas.approximation='none'`.

## Implementation structure

- Reproduction of MOLCAS results with the Ligand Field Theory model
- for a single Dy(III) ion.
- Calculation time: seconds
- Magnetic field
- Single Dy ion
- Real g-tensor
- Rotate the ligand field into the molecular frame
- Liza -this needs more decimal places
- Ligand field parameters (MOLCAS)
- Convert to irreducible spherical tensors
- Supply to Spinach
- Formalism specification

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `dcm2euler()`, `icm2hz()`, `wigner()`, `stev2sph()`, `create()`, `basis()`, `geffect()`, `num2str()`.
