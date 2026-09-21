# kernel/step.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/step.m`
- Signature: `rho=step(spin_system,L,rho,time_step)`
- Total lines: 419

## Purpose

Propagation step function. Computes the action by a matrix exponential without compuing that exponential. Supports one-, two-, and three-point product quadratures. Syntax: rho=step(spin_system,L,rho,time_step)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `comm_series()`, `reordered_taylor()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- L -Liouvillian or Hamiltonian to be used for
- propagation; centre point piecewise-constant
- rule if one matrix is supplied, piecewise-
- linear rule if two matrices {left, right}
- are supplied, piecewise-quadratic if three
- matrices {left, midpoint, right} are given.
- If L is assembled manually from Hamiltonian
- commutation superoperator H, relaxation
- superoperator R, and kinetics superoperator
- K, use L=H+1i*R+1i*K
- State-dependent evolution generators are
- supported: if L{1} is a function handle (see
- iserstep.m documentation), L{2} is current
- time, and L{3} is the method (see iserstep.m
- documentation), the problem is routed to a
- an appropriate Lie group solver.
- rho -state vector or density matrix
- time_step -length of the time step to take

## Outputs

- rho -state vector or density matrix
- Note: we initially had a faithful implementation of the Krylov process
- here -subspace, orthogonalisation, projection, etc., but in all
- our testing it was much inferior to the reordered Taylor process
- that is currently implemented below.
- Note: the peculiar sequence of algebraic operations in the code below
- is designed to minimise the memory footprint in large cases.

## Implementation structure

- Propagation step function. Computes the action by a matrix exponential
- without compuing that exponential. Supports one-, two-, and three-point
- product quadratures. Syntax:
- rho=step(spin_system,L,rho,time_step)
- L - Liouvillian or Hamiltonian to be used for
- propagation; centre point piecewise-constant
- rule if one matrix is supplied, piecewise-
- linear rule if two matrices {left, right}
- are supplied, piecewise-quadratic if three
- matrices {left, midpoint, right} are given.
- If L is assembled manually from Hamiltonian
- commutation superoperator H, relaxation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `iscell()`, `iserstep()`, `grumble()`, `ismember()`, `gpuArray()`, `isergen()`, `cheap_norm()`, `report()`, `num2str()`, `evolution()`, `true()`, `false()`, `comm_series()`, `issparse()`, `cellfun()`, `rho()`.
