# examples/dnp_sol/novel_dnp/novel_field_profile.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/novel_dnp/novel_field_profile.m`
- Signature: `novel_field_profile()`
- Total lines: 83

## Purpose

Field profile of a NOVEL DNP experiment. <I_z> on 1H after a 0.25 us contact time is calculated as a function of electron pulse amplitude and offset. Further information in: Calculation time: seconds.

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: X-band magnet; implemented by `sys.magnet=0.34`.
- Lines 17-18: Electron and two protons; implemented by `sys.isotopes={'E','1H','1H'}`.
- Lines 20-21: Zeeman interactions (g-tensor for trityl, ppm guess for 1H); implemented by `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5],[0 5 0]}`.
- Lines 24-25: Cartesian coordinates; implemented by `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 29-30: Spin temperature; implemented by `inter.temperature=80`.
- Lines 32-33: Hush the output; implemented by `sys.output='hush'`.
- Lines 35-36: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 39-40: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 43-44: Detection state; implemented by `parameters.coil=state(spin_system,'Lz','1H')`.
- Lines 46-47: Experiment parameters; implemented by `parameters.spins={'E','1H'}`.
- Lines 56-57: Microwave offsets, Hz; implemented by `offsets=linspace(-35e6,35e6,71)`.
- Lines 60-61: Loop over offsets; implemented by `parfor m=1:numel(offsets)`.
- Lines 63-64: Localise parameter array; implemented by `localpar=parameters`.
- Lines 66-67: Set the offset; implemented by `localpar.offset=[offsets(m)+reference_point 0]`.
- Lines 69-70: Run the simulation; implemented by `contactcurve=powder(spin_system,@noveldnp,localpar,'esr')`.
- Lines 72-73: Extract the answer; implemented by `field_prof(m)=real(contactcurve(end))`.
- Lines 77-78: Plotting; implemented by `kfigure(); plot(offsets/1e6,field_prof); axis tight`.

### Control flow inferred from the code

- Line 61: `parfor` loop over `m=1:numel(offsets)`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=0.34`.
- Lines 18: computes `sys.isotopes` using `sys.isotopes={'E','1H','1H'}`.
- Lines 21: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5],[0 5 0]}`.
- Lines 22: computes `inter.zeeman.euler` using `inter.zeeman.euler=(pi/180)*{[0 10 0],[0 0 10],[100 0 0]}`.
- Lines 25: computes `inter.coordinates` using `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 30: computes `inter.temperature` using `inter.temperature=80`.
- Lines 33: computes `sys.output` using `sys.output='hush'`.
- Lines 36: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 37: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 40: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 44: computes `parameters.coil` using `parameters.coil=state(spin_system,'Lz','1H')`.
- Lines 47: computes `parameters.spins` using `parameters.spins={'E','1H'}`.
- Lines 48: computes `parameters.timestep` using `parameters.timestep=1e-9`.
- Lines 49: computes `parameters.nsteps` using `parameters.nsteps=250`.
- Lines 50: computes `parameters.flippulse` using `parameters.flippulse=1`.
- Lines 51: computes `parameters.irr_powers` using `parameters.irr_powers=14.48e6`.
- Lines 52: computes `parameters.pulse_dur` using `parameters.pulse_dur=1/(4*parameters.irr_powers)`.
- Lines 53: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.

## Implementation structure

- Field profile of a NOVEL DNP experiment. <I_z> on 1H after a
- 0.25 us contact time is calculated as a function of electron
- pulse amplitude and offset. Further information in:
- Calculation time: seconds.
- X-band magnet
- Electron and two protons
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Cartesian coordinates
- Spin temperature
- Hush the output
- Basis set
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `offsets()`, `powder()`, `field_prof()`, `contactcurve()`, `kfigure()`, `kylabel()`, `kxlabel()`.
