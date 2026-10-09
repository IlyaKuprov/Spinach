# Direct-sum Liouville action on Hilbert matrices

`deriv=hilb_action(spin_system,L,t,rho)` vectorises each physical diagonal substance block, applies the fixed derivative map `L` or the map returned by `L(t,eta)`, and restores a block-diagonal derivative matrix. The Liouville dimension is `sum(D_n^2)`, not the square of the total Hilbert dimension. No inter-species coherence coordinates are allocated.

The formalism must be `zeeman-hilb`, time must be a finite real scalar, and the input density matrix must be finite with the compiled Hilbert dimensions and zero cross-substance entries. Fixed numeric maps must be finite and have the direct-sum Liouville dimensions. The convention is a derivative map without an implicit `-1i`; kinetics and Hilbert IME use this shared operation. Concentrations are not divided out or renormalised.
