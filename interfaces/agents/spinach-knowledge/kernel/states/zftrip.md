# kernel/states/zftrip.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/states/zftrip.m`
- Signature: `rho=zftrip(spin_system,ZFS,pops,Z,B,idx)`
- Total lines: 122

## Purpose

Projection of the zero-field triplet state with user-specified populations of Cartesian ZFS eigenstates onto the higher-field ZFS + Zeeman eigenstates. This is commonly seen in triplet DNP with photo-generated two-electron triplets. Syntax: rho=zftrip(spin_system,ZFS,pops,Z,B,idx)

## Physical / mathematical content

- State-construction utilities. These routines build equilibrium states, singlets, triplets, partner-state expansions, and physically meaningful density operators in the active basis.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- ZFS -3x3 ZFS tensor (Hz) in the laboratory frame
- of reference; use zfs2mat() to get it from
- D, E, and molecular Euler angles
- pops -a three-element vector with populations of
- X, Y, and Z eigenstates of the ZFS tensor
- at zero magnetic field, order: [pX pY pZ]
- X, Y, and Z are labelled using the organic triplet convention |Dzz|>|Dxx|>|Dyy|, under which D and E have opposite signs, -1/3<E/D<0 (Poole, Farach, Jackson, J. Chem. Phys. 61, 2220 (1974), DOI 10.1063/1.1682294); populations quoted in the transition metal convention |Dzz|>|Dyy|>|Dxx| with 0<E/D<1/3 must have X and Y swapped before the call; at E=0 the X and Y states are degenerate (all three when D=0 as well) and the labelling within the degenerate set is undefined; at E/D=-1/3 the X and Z energies are opposite in sign, not degenerate, but equal in magnitude, so the sort by |energy| cannot tell X from Z; in both cases the affected populations must be equal for the result to be meaningful
- Z -3x3 Zeeman interaction tensor (Hz/Tesla) in
- the laboratory frame of reference; use func-
- tions like axrh2mat() to get it from eigen-
- values and molecular Euler angles
- B -magnetic field directed along the Z axis of
- the laboratory frame of reference, Tesla
- idx -index of the electron triplet (use 'E3') in
- the sys.isotopes list

## Outputs

- rho -spin density matrix (Hilbert space) or state
- vector (Liouville space)

## Implementation structure

- Projection of the zero-field triplet state with user-specified
- populations of Cartesian ZFS eigenstates onto the higher-field
- ZFS + Zeeman eigenstates. This is commonly seen in triplet DNP
- with photo-generated two-electron triplets. Syntax:
- rho=zftrip(spin_system,ZFS,pops,Z,B,idx)
- ZFS -3x3 ZFS tensor (Hz) in the laboratory frame
- of reference; use zfs2mat() to get it from
- D, E, and molecular Euler angles
- pops -a three-element vector with populations of
- X, Y, and Z eigenstates of the ZFS tensor
- at zero magnetic field, order: [pX pY pZ]
- X, Y, and Z are labelled using the organic triplet convention |Dzz|>|Dxx|>|Dyy|, under which D and E have opposite signs, -1/3<E/D<0 (Poole, Farach, Jackson, J. Chem. Phys. 61, 2220 (1974), DOI 10.1063/1.1682294); populations quoted in the transition metal convention |Dzz|>|Dyy|>|Dxx| with 0<E/D<1/3 must have X and Y swapped before the call; at E=0 the X and Y states are degenerate (all three when D=0 as well) and the labelling within the degenerate set is undefined; at E/D=-1/3 the X and Z energies are opposite in sign, not degenerate, but equal in magnitude, so the sort by |energy| cannot tell X from Z; in both cases the affected populations must be equal for the result to be meaningful
- Z -3x3 Zeeman interaction tensor (Hz/Tesla) in

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `pauli()`, `ZFS()`, `pops()`, `irr_sph_ten()`, `complex()`, `state()`, `any()`, `isscalar()`, `strcmp()`, `int2str()`.
