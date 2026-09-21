# examples/dnp_sol/steady_state/novel_x_rep_time_ensemble_b1_r.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/novel_x_rep_time_ensemble_b1_r.m`
- Signature: `novel_x_rep_time_ensemble_b1_r()`
- Total lines: 138

## Purpose

Simulation of NOVEL DNP repetition time scan in the steady state with distributions in electron-proton distance and microwave B1 field. Calculation time: hours.

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Simulation of NOVEL DNP repetition time scan in the steady
- state with distributions in electron-proton distance and
- microwave B1 field.
- Calculation time: hours.
- X-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Spin temperature
- Basis set
- Propagator accuracy
- Algorithmic options
- Distance and B1 ensemble

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gaussleg()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `rep_time()`, `dnp_noflip()`, `powder()`, `dnp_flip()`, `kfigure()`, `kylabel()`, `klegend()`, `kxlabel()`, `xlim()`, `savefig()`.
