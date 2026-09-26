# kernel/create.m

- Signature: `spin_system=create(sys,inter)`

## Purpose

The entry function of the Spinach kernel that creates the spin system object that the rest of the library requires to run. It checks and absorbs interaction specifications, and writes some useful diagnostics to the console. Syntax: spin_system=create(sys,inter)

## Physical / mathematical content

- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Parameters / inputs

- sys -spin system and instrument specification
- structure, see the spin system specifica-
- tion section of the online manual
- inter -interaction specification structure, see
- see the spin system specification section
- of the online manual

## Outputs

- spin_system -the primary object used by Spinach
- to store simulation information
- Note: inter.modes.carriers declares, for each bosonic mode, the labo-
- ratory frequency of the rotating frame in which inter.modes.frqs
- is specified; declared frequencies are then detunings and may be
- negative, and thermal occupations are computed from the physical
- frequency, meaning the sum of the carrier and the detuning.
- Note: inter.modes.t2_times is interpreted at the declared temperature:
- the pure dephasing rate is extracted as 1/T2-kappa*(1+2*nbar)/2,
- where kappa is the amplitude damping rate and nbar the thermal
- occupation at the physical mode frequency.
- Note: quadrature operators of bosonic modes follow the (a+a')/sqrt(2)
- normalisation everywhere, in inter.modes.longitudinal as well as
- in the modulation channels.

## Implementation structure

- The entry function of the Spinach kernel that creates the spin
- system object that the rest of the library requires to run. It
- checks and absorbs interaction specifications, and writes some
- useful diagnostics to the console. Syntax:
- spin_system=create(sys,inter)
- sys -spin system and instrument specification
- structure, see the spin system specifica-
- tion section of the online manual
- inter -interaction specification structure, see
- see the spin system specification section
- of the online manual
- spin_system -the primary object used by Spinach
