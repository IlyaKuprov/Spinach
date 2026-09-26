# examples/dnp_sol/solid_effect_timedep_2.m

- Signature: `solid_effect_timedep_2()`

## Purpose

A simulation of solid effect DNP for a tilted linear chain of protons positioned at distances of 4+n^2 Angstrom with n=2:6 and the electron located at the origin. Weizmann DNP relaxation model is used with se- cond order Krylov-Bogolyubov average Hamiltonian theory and state space restriction to five-spin orders. Calculation time: minutes with a Tesla A100 GPU

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- A simulation of solid effect DNP for a tilted linear chain of protons
- positioned at distances of 4+n^2 Angstrom with n=2:6 and the electron
- located at the origin. Weizmann DNP relaxation model is used with se-
- cond order Krylov-Bogolyubov average Hamiltonian theory and state
- space restriction to five-spin orders.
- Calculation time: minutes with a Tesla A100 GPU
- Magnetic field
- Electron
- Nuclei
- Relaxation theory
- Microwave power and offset
- Basis set
