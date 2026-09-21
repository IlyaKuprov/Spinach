# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r.m`
- Signature: `xix_q_field_profile_ensemble_r()`
- Total lines: 96

## Purpose

Simulation of XiX DNP field profile in the steady state with electron-proton distance ensemble averaging. Calculation time: minutes.

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Simulation of XiX DNP field profile in the steady state
- with electron-proton distance ensemble averaging.
- Calculation time: minutes.
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Spin temperature
- Basis set
- Propagator accuracy
- Algorithmic options
- Distance ensemble, Gauss-Legendre points
- Microwave resonance offsets, Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gaussleg()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `dnp()`, `powder()`, `kfigure()`, `kylabel()`, `kxlabel()`, `savefig()`.
