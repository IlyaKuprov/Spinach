# examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r.m`
- Signature: `xix_q_con_time_ensemble_r()`
- Total lines: 111

## Purpose

Simulation of XiX DNP contact time dependence in the steady state with electron-proton distance ensembles. Calculation time: minutes.

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Simulation of XiX DNP contact time dependence in the
- steady state with electron-proton distance ensembles.
- Calculation time: minutes.
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Spin temperature
- Basis set
- Propagator accuracy
- Algorithmic options
- Distance ensemble
- XiX loop count

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gaussleg()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `loop_counts()`, `dnp()`, `powder()`, `kfigure()`, `kylabel()`, `kxlabel()`, `savefig()`.
