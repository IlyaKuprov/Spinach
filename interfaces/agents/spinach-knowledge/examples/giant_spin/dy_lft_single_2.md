# examples/giant_spin/dy_lft_single_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/dy_lft_single_2.m`
- Signature: `dy_lft_single_2()`
- Total lines: 91

## Purpose

A demonstration that most lanthanide complexes are in the ZFS limit for the purposes of relaxation theory. One of the figures from our forthcoming papers on the subject. Calculation time: hours.

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnetic field; implemented by `sys.magnet=1.0`.
- Lines 15-16: Single Dy ion; implemented by `sys.isotopes={'E16'}`.
- Lines 18-19: Real g-tensor; implemented by `D=[ 1.322766699 1.324261429 1.328750739]`.
- Lines 25-26: Rotate the ligand field into the molecular frame; implemented by `R=[ 0.56699461681939 -0.81924972266747 0.08571462189795`.
- Lines 31-32: Ligand field parameters (MOLCAS); implemented by `Bkq{2}=[ 0.14365204369342E-01`.
- Lines 60-61: Convert to irreducible spherical tensors; implemented by `for k=[2 4 6]`.
- Lines 66-69: Supply to Spinach; implemented by `inter.giant.coeff={{[0 0 0],SBkq{2}, [0 0 0 0 0 0 0],SBkq{4}, [0 0 0 0 0 0 0 0 0 0 0],SBkq{6}}}`.
- Lines 73-74: Formalism specification; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 77-78: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 81-82: Experiment parameters; implemented by `parameters.fields=[0 500]`.
- Lines 87-88: Run the simulation; implemented by `fieldscan_enlev(spin_system,parameters)`.

### Control flow inferred from the code

- Line 61: `for` loop over `k=[2 4 6]`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=1.0`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'E16'}`.
- Lines 19: computes `D` using `D=[ 1.322766699 1.324261429 1.328750739]`.
- Lines 20: computes `V` using `V=[ 0.566995 -0.819250 0.085715`.
- Lines 23: computes `inter.zeeman.matrix` using `inter.zeeman.matrix={V'*diag(D)*V}`.
- Lines 26: computes `R` using `R=[ 0.56699461681939 -0.81924972266747 0.08571462189795`.
- Lines 29: computes `[alp,bet,gam]` using `[alp,bet,gam]=dcm2euler(R')`.
- Lines 32: computes `Bkq{2}` using `Bkq{2}=[ 0.14365204369342E-01`.
- Lines 37: computes `Bkq{4}` using `Bkq{4}=[-0.32862425456238E-01`.
- Lines 46: computes `Bkq{6}` using `Bkq{6}=[-0.61447247154955E-04`.
- Lines 62: computes `Bkq{k}` using `Bkq{k}=icm2hz(Bkq{k})`.
- Lines 63: computes `SBkq{k}` using `SBkq{k}=wigner(k,alp,bet,gam)*stev2sph(k,Bkq{k})`.
- Lines 67-69: computes `inter.giant.coeff` using `inter.giant.coeff={{[0 0 0],SBkq{2}, [0 0 0 0 0 0 0],SBkq{4}, [0 0 0 0 0 0 0 0 0 0 0],SBkq{6}}}`.
- Lines 70-71: computes `inter.giant.euler` using `inter.giant.euler={{[0 0 0],[0 0 0],[0 0 0], [0 0 0],[0 0 0],[0 0 0]}}`.
- Lines 74: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 75: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 78: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 82: computes `parameters.fields` using `parameters.fields=[0 500]`.

## Implementation structure

- A demonstration that most lanthanide complexes are in the
- ZFS limit for the purposes of relaxation theory. One of the
- figures from our forthcoming papers on the subject.
- Calculation time: hours.
- Magnetic field
- Single Dy ion
- Real g-tensor
- Rotate the ligand field into the molecular frame
- Ligand field parameters (MOLCAS)
- Convert to irreducible spherical tensors
- Supply to Spinach
- Formalism specification

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `dcm2euler()`, `icm2hz()`, `wigner()`, `stev2sph()`, `create()`, `basis()`, `fieldscan_enlev()`.
