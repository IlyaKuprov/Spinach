# examples/dnp_sol/tppm_dnp/tppm_contact_curve.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/tppm_dnp/tppm_contact_curve.m`
- Signature: `tppm_contact_curve()`
- Total lines: 61

## Purpose

The transformation of -E_z into I_z during the contact time of the TPPM DNP experiment. Further information in: Calculation time: seconds

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Q-band magnet; implemented by `sys.magnet=1.2142`.
- Lines 16-17: Electron and two protons; implemented by `sys.isotopes={'E','1H','1H'}`.
- Lines 19-20: Zeeman interactions (g-tensor for trityl, ppm guess for 1H); implemented by `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5],[0 5 0]}`.
- Lines 23-24: Cartesian coordinates; implemented by `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 28-29: Spin temperature; implemented by `inter.temperature=80`.
- Lines 31-32: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Detection state; implemented by `parameters.coil=state(spin_system,'Lz','1H')`.
- Lines 42-43: Experiment parameters; implemented by `parameters.spins={'E','1H'}`.
- Lines 52-53: Run the calculation; implemented by `contact_curve=powder(spin_system,@xixdnp,parameters,'esr')`.
- Lines 55-56: Plotting; implemented by `time_axis=linspace(0,2*parameters.pulse_dur*parameters.nloops,parameters.nloops+1)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=1.2142`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'E','1H','1H'}`.
- Lines 20: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5],[0 5 0]}`.
- Lines 21: computes `inter.zeeman.euler` using `inter.zeeman.euler=(pi/180)*{[0 10 0],[0 0 10],[100 0 0]}`.
- Lines 24: computes `inter.coordinates` using `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 29: computes `inter.temperature` using `inter.temperature=80`.
- Lines 32: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 33: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `parameters.coil` using `parameters.coil=state(spin_system,'Lz','1H')`.
- Lines 43: computes `parameters.spins` using `parameters.spins={'E','1H'}`.
- Lines 44: computes `parameters.offset` using `parameters.offset=[(-13+4)*1e6 0]`.
- Lines 45: computes `parameters.irr_powers` using `parameters.irr_powers=33e6`.
- Lines 46: computes `parameters.pulse_dur` using `parameters.pulse_dur=16e-9`.
- Lines 47: computes `parameters.nloops` using `parameters.nloops=250`.
- Lines 48: computes `parameters.phase` using `parameters.phase=2*pi/3`.
- Lines 49: computes `parameters.grid` using `parameters.grid='rep_2ang_400pts_sph'`.
- Lines 50: computes `parameters.needs` using `parameters.needs={'aniso_eq'}`.

## Implementation structure

- The transformation of -E_z into I_z during the contact time
- of the TPPM DNP experiment. Further information in:
- Calculation time: seconds
- Q-band magnet
- Electron and two protons
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Cartesian coordinates
- Spin temperature
- Basis set
- Spinach housekeeping
- Detection state
- Experiment parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `powder()`, `kfigure()`, `kxlabel()`, `kylabel()`.
