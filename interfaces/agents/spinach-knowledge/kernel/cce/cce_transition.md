# kernel/cce/cce_transition.m

## Native role and invocation

Explicit observed-transition selection utility shared with native model compilation. Selection does not itself prepare the observed state.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `[transition,meta]=cce_transition(spin_system,parameters,assumptions)`. Engine records are not public experiment parameter structures.

[MATLAB source](../../../../../kernel/cce/cce_transition.m)

## Contract

observed and transition.spins are identical ordered physical spin-ID rows. Select exactly one of ms (n-observed-by-2 real magnetic quantum numbers), levels (two distinct ascending-energy indices), or orthonormal alpha/beta columns. Product states use descending magnetic-quantum-number ordering; the first listed spin is the first tensor factor. Explicit states specify the physical observed transition.

The caller's global random-stream state is preserved, including on error. Only the observed subsystem is retained and given an exact zeeman-hilb basis. Stock assume, hamiltonian, and orientation are used at the actual field and the supplied radian Euler-angle row, or the reported zero-angle frame. Canonical tensors already use rad/s. Energy-index selectors reject either selected state's unresolved degeneracy at the dimension-scaled floating-point resolution; explicitly chosen states remain possible. For nondegenerate eigenstates the largest component has positive real phase.

transition returns spins, alpha, and beta. meta gives the selector, ordered multiplicities, frame, assumptions, expectation energies, signed beta-minus-alpha gap, nearest gaps, residual norms, and degeneracy threshold, all frequency quantities in rad/s. ms selectors need not be exact eigenstates; residuals expose that distinction, nearest_gaps is NaN for non-eigenstates, and eigenvalues reports the observed spectrum. No state population is prepared, and no secular projection is silently applied.
The first argument may also be a compiled `cce_model`. For product-ms selectors in the coherent mixed-product regime without bath jumps this reuses an observed-only native retained generator with mean field disabled. Canonical create inputs, Liouville representations, sampled states, and dissipative systems retain the stock construction. Ordering, energies, degeneracy handling, and residual metadata follow the same selector contract.

The entry validator accepts either a canonical Spinach system or the compiled-model form and checks the same transition selection before observed-subsystem construction.
