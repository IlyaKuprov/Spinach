# examples/dnp_sol/xix_dnp/xix_parameter_scan.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/xix_dnp/xix_parameter_scan.m`
- Signature: `xix_parameter_scan()`
- Total lines: 100

## Purpose

2D parameter scan of a XiX DNP experiment. <I_z> after a set contact time is calculated as a function of electron pulse amplitude and offset. Further information in: Calculation time: minutes (a large powder grid is needed).

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- 2D parameter scan of a XiX DNP experiment. <I_z> after a set
- contact time is calculated as a function of electron pulse
- amplitude and offset. Further information in:
- Calculation time: minutes (a large powder grid is needed).
- Q-band magnet
- Electron and two protons
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Cartesian coordinates
- Spin temperature
- Basis set
- Spinach housekeeping
- Detection state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `kfigure()`, `offsets()`, `nutfrqs()`, `powder()`, `dnp_surf()`, `contact_curve()`, `contourf()`, `kylabel()`, `kxlabel()`, `kcolourbar()`.
