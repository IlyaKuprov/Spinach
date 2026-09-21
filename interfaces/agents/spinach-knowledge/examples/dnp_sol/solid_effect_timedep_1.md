# examples/dnp_sol/solid_effect_timedep_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/solid_effect_timedep_1.m`
- Signature: `solid_effect_timedep_1()`
- Total lines: 91

## Purpose

A simulation of solid effect DNP for a tilted linear chain of three protons positioned at distances 7, 10 and 14 Angstrom from electron located at the origin. Weizmann DNP relaxation model is used with second order Krylov-Bogolyubov average Hamiltonian theory. Calculation time: seconds

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- A simulation of solid effect DNP for a tilted linear chain of three
- protons positioned at distances 7, 10 and 14 Angstrom from electron
- located at the origin. Weizmann DNP relaxation model is used with
- second order Krylov-Bogolyubov average Hamiltonian theory.
- Calculation time: seconds
- Magnetic field
- Spin system
- Relaxation theory
- Microwave power and offset
- Basis set
- Spinach housekeeping
- Experiment parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `euler2dcm()`, `create()`, `basis()`, `solid_effect()`, `kfigure()`, `scale_figure()`, `subplot()`, `answer()`, `kylabel()`, `kxlabel()`, `klegend()`, `set()`.
