# kernel/states/equilibrium.m

- MATLAB implementation: [kernel/states/equilibrium.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/states/equilibrium.m)

## Purpose

Return the Boltzmann thermal-equilibrium density matrix (Hilbert space) or state vector (Liouville space) for the configured spin-system temperature. The full interface is `rho=equilibrium(spin_system,I,Q,euler_angles)`, but the accepted call forms are one, two, or four arguments:

- `rho=equilibrium(spin_system)` builds the Hamiltonian internally under the `'labframe'` assumption and uses its isotropic part.
- `rho=equilibrium(spin_system,I)` uses the supplied isotropic Hamiltonian only.
- `rho=equilibrium(spin_system,I,Q,euler_angles)` adds the anisotropic contribution at the requested orientation.

A three-argument call is not implemented.

## Inputs and constraints

- Multi-substance Zeeman systems use local vectorised identities in Liouville space and independently trace-normalised blocks in Hilbert space. All three call forms retain the single-substance algorithms.
- `spin_system` must have a temperature in `spin_system.rlx.temperature` (configured through `inter.temperature`); an empty or exactly zero temperature is rejected. The routine uses `hbar/(kbol*T)` from the configured physical constants.
- `I` is numeric: the isotropic Hamiltonian in Hilbert space, or its **left-side product superoperator** in Liouville space. In Liouville calculations, neither `I` nor anisotropic terms may be supplied as commutation superoperators.
- For the four-argument form, `Q` is a cell array of anisotropic Hamiltonian terms and `euler_angles` is a real three-element vector in radians, giving the system orientation relative to the input orientation. The orientation-dependent term is added to `I`.
- In segmented systems, the assembled Hamiltonian (after adding the orientation-dependent contribution, when supplied) must be block diagonal by substance. Any nonzero cross-substance block raises `Spinach:equilibrium:crossSubstanceHamiltonian` before propagation.
- Hamiltonians `I` and `Q` generated with `hamiltonian.m` must use the `'labframe'` assumption.

## Numerical mechanism and limitations

The calculation forms the imaginary-time thermal state with `beta=hbar/(kbol*T)`. Liouville formalisms propagate the thermodynamic unit state with `step` at `-1i*beta` and normalise each substance block by its own overlap with that unit state; each normalised block is then multiplied by its concentration. The thermal unit used for imaginary-time propagation is unweighted, so zero-concentration substances need no division by zero. Zeeman Hilbert formalism obtains an imaginary-time propagator, normalises by the trace, and uses scaling-and-squaring with intermediate normalisation/cleanup for numerical stability, then multiplies the final density matrix by the concentration. For a Hilbert direct sum, each diagonal block runs that same calculation independently before block-diagonal assembly. Before propagation, the action of each diagonal Hamiltonian block on its own unit state is checked independently with the existing `1e-10` one-norm threshold. Exactly zero Hamiltonian blocks in segmented Liouville systems are valid left products and are exempt from this check, whether spinful or spin-free; their equilibrium is the concentration-weighted local unit state. A vanishing action in a checked block raises `Spinach:equilibrium:notLeftProduct` and identifies the substance; a valid left-product block cannot hide a commutator in another block. This is a rejection test for vanishing unit action, not a complete algebraic certification of an arbitrary supplied matrix. In Liouville space, a NaN result is reported as a too-low-temperature accuracy failure, for which the source error message directs the caller to switch to Hilbert space. Ground-state degeneracy is common, so absolute-zero equilibrium is unsupported.

Source and credits: [Spinach Wiki: equilibrium.m](https://spindynamics.org/wiki/index.php?title=equilibrium.m). Luke Edwards and Ilya Kuprov (contact details are in the source file).

Wavefunction equilibrium is explicitly rejected by `Spinach:equilibrium:wavefunction`; a finite-temperature mixed state is not a single ket.
