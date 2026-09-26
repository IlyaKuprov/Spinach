# interfaces/g2spinach.m

- Signature: `[sys,inter]=g2spinach(props,particles,references,options)`

## Purpose

Makes Spinach data structures from parsed outputs of electronic structure theory packages, such as Gaussian and ORCA. Syntax: [sys,inter]=g2spinach(props,particles,references,options)

## Physical / mathematical content

- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

## Parameters / inputs

- props -the output of gparse() function
- particles -a cell array of the following form:
- {{'H','1H'},{'N','15N'}...}
- giving the list of elements and isotopes that
- should be imported. If the isotope list contains
- an electron, e.g. {{'E','E'},{'H','1H'}...},
- then EPR mode is assumed -chemical shielding
- and scalar couplings are ignored, but g-tensor
- and hyperfine couplings are included.
- references -a vector of absolute shielding values for
- the reference substances that are to be
- placed at zero ppm chemical shift; you ne-
- ed to run separate electronic structure
- theory calculations for those substances
- wuth the same method. Absolute isotropic
- shielding values for tetramethylsilane in
- vacuum are:
- GIAO 13C 1H
- B3LYP/6-31G* 189.6621 32.1833
- B3LYP/6-311+G(2d,p) 182.4485 31.8201
- HF/6-31G* 199.9711 32.5957
- HF/6-311+G(2d,p) 192.5828 32.0710
- CSGT 13C 1H
- B3LYP/6-31G* 188.5603 29.1952
- B3LYP/6-311+G(2d,p) 182.1386 31.7788
- HF/6-31G* 196.8670 29.5517
- HF/6-311+G(2d,p) 192.5701 31.5989
- This setting is ignored when electrons are
- present in the system.
- options.min_j -scalar coupling threshold in Hz. J-coup-
- lings smaller than this value will be
- ignored in the NMR mode.
- options.min_hfc -hyperfine coupling threshold in Hz. Hy-
- perfine tensors with a Frobenius norm
- smaller than this value will be ignored
- in the EPR mode.
- options.purge -if set to 'on' in EPR mode, removes the
- spins with hyperfine coupling below
- options.min_hfc from the spin system.
- options.no_xyz -if set to 1, causes the function to ig-
- nore the coordinate information and
- only keep the interaction tensors

## Outputs

- sys.isotopes Nspins x 1 cell array of strings
- inter.coordinates Nspins x 3 dense matrix, Angstrom. Not
- returned if there is an electron in the
- isotope list (in the EPR case it is not
- a good idea to use the molecular
- coordinates for spins).
- inter.zeeman.matrix Nspins x 1 cell array of 3x3 matrices,
- ppm for nuclei, g-tensor for electrons.
- Zero interactions have zero matrices.
- inter.coupling.matrix Nspins x Nspins cell array of 3x3 mat-
- rices, all in Hz. Zero interactions have
- zero matrices.
- inter.coupling.scalar Nspins x Nspins cell array of scalar
- couplings, all in Hz. Zero couplings are
- returned as zeros.
- inter.spinrot.matrix spin-rotation coupling tensors for
- each nucleus

## Implementation structure

- Makes Spinach data structures from parsed outputs of electronic
- structure theory packages, such as Gaussian and ORCA. Syntax:
- [sys,inter]=g2spinach(props,particles,references,options)
- props -the output of gparse() function
- particles -a cell array of the following form:
- {{'H','1H'},{'N','15N'}...}
- giving the list of elements and isotopes that
- should be imported. If the isotope list contains
- an electron, e.g. {{'E','E'},{'H','1H'}...},
- then EPR mode is assumed -chemical shielding
- and scalar couplings are ignored, but g-tensor
- and hyperfine couplings are included.
