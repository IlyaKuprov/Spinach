# kernel/average.m

- Signature: `H=average(spin_system,Hp,H0,Hm,omega,theory)`

## Purpose

Average Hamiltonian theories under Zeeman interaction rotating frame transformations. Syntax: H=average(spin_system,Hp,H0,Hm,omega,theory)

## Physical / mathematical content

- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Parameters / inputs

- Hp -the part of the rotating frame Hamiltonian that has positive
- frequency +omega under the rotating frame transformation
- H0 -the part of the rotating frame Hamiltonian that has zero
- frequency under the rotating frame transformation
- Hm -the part of the rotating frame Hamiltonian that has negative
- frequency -omega under the rotating frame transformation
- omega -the frequency of the rotating frame transformation, rad/s
- theory -the level of the average Hamiltonian theory:
- 'ah_first_order' -first order in Waugh theory
- 'ah_second_order' -second order in Waugh theory
- 'ah_third_order' -third order in Waugh theory
- 'matrix_log' -exact algorithm (very expensive,
- uses dense matrix algebra)
- 'kb_first_order' -first order in Krylov-Bogolyubov
- theory (DNP experiments only)
- 'kb_second_order' -second order in Krylov-Bogolyubov
- theory (DNP experiments only)
- 'kb_third_order' -third order in Krylov-Bogolyubov
- theory (DNP experiments only)

## Outputs

- H -average Hamiltonian
- Note: Krylov-Bogolyubov averging theory as applied to DNP systems is
- described in detail here:

## Implementation structure

- Average Hamiltonian theories under Zeeman interaction rotating frame
- transformations. Syntax:
- H=average(spin_system,Hp,H0,Hm,omega,theory)
- Hp - the part of the rotating frame Hamiltonian that has positive
- frequency +omega under the rotating frame transformation
- H0 - the part of the rotating frame Hamiltonian that has zero
- frequency under the rotating frame transformation
- Hm - the part of the rotating frame Hamiltonian that has negative
- frequency -omega under the rotating frame transformation
- omega - the frequency of the rotating frame transformation, rad/s
- theory - the level of the average Hamiltonian theory:
- 'ah_first_order' -first order in Waugh theory
