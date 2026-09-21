# examples/dnp_sol/steady_state/tppm_q_con_time_ensemble_b1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/tppm_q_con_time_ensemble_b1.m`
- Signature: `tppm_q_con_time_ensemble_b1()`
- Total lines: 117

## Purpose

Simulation of TPPM DNP contact time dependence in the steady state with electron Rabi frequency ensemble. Calculation time: hours.

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Simulation of TPPM DNP contact time dependence in the
- steady state with electron Rabi frequency ensemble.
- Calculation time: hours.
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Spin temperature
- Cartesian coordinates
- Get electron-nuclear distance
- Basis set
- Propagator accuracy
- Algorithmic options

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cell2mat()`, `xyz()`, `gaussleg()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `loop_counts()`, `dnp()`, `powder()`, `kfigure()`, `kylabel()`, `kxlabel()`, `savefig()`.
