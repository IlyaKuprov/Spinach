# kernel/kinetics/chem_concs.m

`concs=chem_concs(spin_system,eta)` reads the trace functional of each substance, using the compiled offsets: the spherical-tensor unit coordinate, the vectorised identity in Zeeman Liouville space, or the local matrix trace in Hilbert space. In either Liouville formalism, the finite column `eta` contains one or more complete spin blocks in space-times-spin order; the result is an `nvoxels`-by-`nsubst` array. In Hilbert space, `eta` is a finite block-diagonal matrix and the result is a one-row concentration array. Nonzero inter-species coherences are rejected, not discarded.

Concentrations are state coordinates, including zero-population and spin-free blocks. The function neither clips nor renormalises them. Wavefunctions are explicitly rejected because a pure ket has no mixed-state concentration trace functional.
