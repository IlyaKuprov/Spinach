# examples/spin_chemistry/cidnp_geminate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/cidnp_geminate.m`
- Signature: `cidnp_geminate()`
- Total lines: 59

## Purpose

A basic example of the geminate CIDNP effect simulation. Calculation time: seconds

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- A basic example of the geminate CIDNP effect simulation.
- Calculation time: seconds
- System specification
- Basis set
- Spinach housekeeping
- Get the Hamiltonian
- Get the kinetics superoperator
- Get the initial state
- Double up the problem (no dynamics assumed in the product subspace)
- Set up a reaction ("whatever is leaving reactants must appear in products")
- Assemble the Liouvillian
- Evolve for a microsecond

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `kinetics()`, `singlet()`, `evolution()`, `state()`, `rho()`, `num2str()`.
