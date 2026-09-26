# kernel/states/equilibrium.m

- Signature: `rho=equilibrium(spin_system,I,Q,euler_angles)`

## Purpose

Returns the thermal equilibrium state at the current temperature. If the anisotropic part and the orientation parameters are not given, uses the isotropic Hamiltonian, otherwise uses the full Hamiltonian at the speci- fied orientation. Syntax: rho=equilibrium(spin_system,I,Q,euler_angles)

## Physical / mathematical content

- State-construction utilities. These routines build equilibrium states, singlets, triplets, partner-state expansions, and physically meaningful density operators in the active basis.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- I -isotropic part of the Hamiltonian left side pro-
- duct superoperator (in Liouville space) or Hamil-
- tonian (in Hilbert space). If this argument is
- omitted, the Hamiltonian is built here and used
- to compute the thermal equilibrium state.
- Q -irreducible components of the anisotropic part
- of the Hamiltonian left side product superopera-
- tor (in Liouville space) or Hamiltonian (in Hil-
- bert space), as returned by hamiltonian.m; this
- is needed when the thermal equilibrium state de-
- pends on the system orientation.
- euler_angles -a row vector of Euler angles (in radians) speci-
- fying the system orientation relative to the in-
- put orientation. If the angles are not supplied,
- only isotropic part of the Hamiltonian is used.

## Outputs

- rho -thermal equilibrium density matrix (Hilbert spa-
- ce) or state vector (Liouville space).
- WARNING: Liouville space calculations must supply left side product su-
- peroperators, not commutation superoperators.
- WARNING: assumptions supplied to the hamiltonian.m call that generates
- I and Q must be 'labframe'.
- WARNING: spin system ground states are commonly degenerate; absolute
- zero temperatures are not supported.

## Implementation structure

- Returns the thermal equilibrium state at the current temperature. If the
- anisotropic part and the orientation parameters are not given, uses the
- isotropic Hamiltonian, otherwise uses the full Hamiltonian at the speci-
- fied orientation. Syntax:
- rho=equilibrium(spin_system,I,Q,euler_angles)
- I -isotropic part of the Hamiltonian left side pro-
- duct superoperator (in Liouville space) or Hamil-
- tonian (in Hilbert space). If this argument is
- omitted, the Hamiltonian is built here and used
- to compute the thermal equilibrium state.
- Q -irreducible components of the anisotropic part
- of the Hamiltonian left side product superopera-
