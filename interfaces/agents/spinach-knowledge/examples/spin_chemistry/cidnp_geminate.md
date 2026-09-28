# examples/spin_chemistry/cidnp_geminate.m

- Signature: `cidnp_geminate()`

## Purpose

A basic example of the geminate CIDNP effect simulation. Calculation time: seconds

## Physical / mathematical content

- At 14.1 T the system contains two electron spins and one 1H spin, with g factors `2.0023`, `2.0024`, and `1.0`. The electron pair has a singlet initial state and couples to the proton with `J=1e7`; recombination uses Haberkorn theory with rates `[1e7 0]`.
- The model duplicates the state space into reactant and product sectors; the reaction superoperator removes population from reactants and transfers it to products, where no further dynamics are assumed.

## Numerical / algorithmic content

- It evolves the assembled Liouvillian for `1e-6 s` and reports the proton `Lz` magnetisation in the reactant and product sectors separately.

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
