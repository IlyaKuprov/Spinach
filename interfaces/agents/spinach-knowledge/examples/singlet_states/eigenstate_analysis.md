# examples/singlet_states/eigenstate_analysis.m

- Signature: `eigenstate_analysis()`

## Purpose

Analyses the part of an allyl-pyruvate singlet state that is stationary under the drift Hamiltonian, then repeats the projection after applying a proton offset and spin-lock term. It is a Hamiltonian/eigenstate analysis, not a relaxation or lifetime calculation.

## Spin system and Hamiltonian

The source obtains a fitted `1H/13C` system from `allyl_pyruvate({'1H','13C'})`, sets the field to `14.1 T`, dilutes to the `13C` isotopomer, and selects subsystem 4. It uses the singlet on spin labels 3 and 4, an unapproximated `zeeman-hilb` basis, and the isotropic NMR Hamiltonian. The full Hamiltonian is explicitly symmetrised before diagonalisation.

## Stationary-state projections

For the singlet operator, the code reports its norm, applies `remncomm` using the Hamiltonian eigenvectors and eigenvalues to retain the commuting component, removes the identity component with `remtrace`, and reports that norm. It then removes the normalised `Lz–Lz` component on spins 3 and 4 and reports the remaining norm. These are reported projections; the source contains no printed numerical results here.

The second pass sets the proton offset parameter to `2850` (the source does not state its unit) and adds the documented `1 kHz` proton spin-lock term, `2*pi*1000*Lx`, before repeating the projections. No gradients, relaxation, storage interval, or time evolution are specified.

## Source

[examples/singlet_states/eigenstate_analysis.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/singlet_states/eigenstate_analysis.m)
