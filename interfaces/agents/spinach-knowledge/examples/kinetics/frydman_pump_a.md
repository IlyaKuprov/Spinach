# examples/kinetics/frydman_pump_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/kinetics/frydman_pump_a.m`
- Signature: `frydman_pump_a()`
- Total lines: 132

## Purpose

Lucio Frydman's water exchange based spin-lock pump, Figure 2 from https://doi.org/10.1016/j.jmr.2021.107083 Calculation time: seconds.

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file also defines local helper function(s): `frydman_pump()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Implementation structure

- Lucio Frydman's water exchange based spin-lock pump, Figure 2
- from https://doi.org/10.1016/j.jmr.2021.107083
- Calculation time: seconds.
- Number of water protons
- Magnet field
- Core spin system
- Add water protons
- Chemical shifts
- Scalar couplings, Hz
- Estimated relaxation times, seconds
- Relaxation theory
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `liquid()`, `state()`, `kfigure()`, `scale_figure()`, `subplot()`, `klegend()`, `kylabel()`, `kxlabel()`, `frydman_pump()`, `operator()`, `equilibrium()`, `dictum()`, `hamiltonian()`.
