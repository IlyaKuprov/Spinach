# examples/dnp_sol/steady_state/xix_w_pulse_dur_single.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/xix_w_pulse_dur_single.m`
- Signature: `xix_w_pulse_dur_single()`
- Total lines: 102

## Purpose

2D parameter scan of XiX DNP in the steady state. Calculation time: hours.

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- 2D parameter scan of XiX DNP in the steady state.
- Calculation time: hours.
- W-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Spin temperature
- Cartesian coordinates
- Get electron-nuclear distance
- Basis set
- Propagator accuracy
- Algorithmic options
- Electron pulse duration grid, s

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cell2mat()`, `xyz()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `kfigure()`, `pulse_durs()`, `dnp()`, `powder()`, `set()`, `kylabel()`, `kcolourbar()`, `kxlabel()`, `savefig()`.
