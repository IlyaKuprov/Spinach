# kernel/kinetics/chem_concs.m

`concs=chem_concs(spin_system,eta)` reads the spherical-tensor unit coordinate of each substance, using the compiled offsets. The finite column `eta` contains one or more complete spin blocks in space-times-spin order; the result is an `nvoxels`-by-`nsubst` array.

Concentrations are state coordinates, including zero-population and spin-free blocks. The function neither clips nor renormalises them. Other formalisms are explicitly rejected until their trace-functional implementation is supplied.
