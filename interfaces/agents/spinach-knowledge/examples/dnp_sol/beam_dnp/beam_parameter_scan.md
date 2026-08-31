# examples/dnp_sol/beam_dnp/beam_parameter_scan.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/beam_dnp/beam_parameter_scan.m`
- Signature: `beam_parameter_scan()`
- Total lines: 99

## Purpose

2D parameter scan of a BEAM DNP experiment. <I_z> after a set contact time is calculated as a function of electron pulse amplitude and offset. Further information in: Calculation time: minutes (a large powder grid is needed).

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: X-band magnet; implemented by `sys.magnet=0.3483`.
- Lines 17-18: Electron and two protons; implemented by `sys.isotopes={'E','1H','1H'}`.
- Lines 20-21: Zeeman interactions (g-tensor for trityl, ppm guess for 1H); implemented by `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5],[0 5 0]}`.
- Lines 24-25: Cartesian coordinates; implemented by `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 29-30: Spin temperature; implemented by `inter.temperature=80`.
- Lines 32-33: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Detection state; implemented by `parameters.coil=state(spin_system,'Lz','1H')`.
- Lines 43-44: Experiment parameters; implemented by `parameters.spins={'E','1H'}`.
- Lines 50-51: MW resonance offset grid, Hz; implemented by `offsets=linspace(-60e6,60e6,120)`.
- Lines 54-55: Electron pulse nutation frequency grid, Hz; implemented by `nutfrqs=linspace(20e6,40e6,30)`.
- Lines 57-58: Hush the output; implemented by `spin_system.sys.output='hush'`.
- Lines 60-61: Generate axis ticks; implemented by `[X,Y]=meshgrid(offsets/1e6,nutfrqs/1e6)`.
- Lines 63-64: Get the figure going; implemented by `kfigure(); dnp_surf=zeros(numel(nutfrqs),numel(offsets))`.
- Lines 66-67: Loop over offsets; implemented by `for m=1:numel(offsets)`.
- Lines 69-70: Set the offsets; implemented by `parameters.offset=[offsets(m)+reference_point 0]`.
- Lines 72-73: Parallel loop over nutation freqs; implemented by `parfor n=1:numel(nutfrqs)`.
- Lines 75-76: Localise parameter array; implemented by `localpar=parameters`.

### Control flow inferred from the code

- Line 67: `for` loop over `m=1:numel(offsets)`.
- Line 73: `parfor` loop over `n=1:numel(nutfrqs)`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=0.3483`.
- Lines 18: computes `sys.isotopes` using `sys.isotopes={'E','1H','1H'}`.
- Lines 21: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5],[0 5 0]}`.
- Lines 22: computes `inter.zeeman.euler` using `inter.zeeman.euler=(pi/180)*{[0 10 0],[0 0 10],[100 0 0]}`.
- Lines 25: computes `inter.coordinates` using `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 30: computes `inter.temperature` using `inter.temperature=80`.
- Lines 33: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `parameters.coil` using `parameters.coil=state(spin_system,'Lz','1H')`.
- Lines 44: computes `parameters.spins` using `parameters.spins={'E','1H'}`.
- Lines 45: computes `parameters.pulse_dur` using `parameters.pulse_dur=[20.0e-9 28.7e-9]`.
- Lines 46: computes `parameters.nloops` using `parameters.nloops=165`.
- Lines 47: computes `parameters.grid` using `parameters.grid='rep_2ang_800pts_sph'`.
- Lines 48: computes `parameters.needs` using `parameters.needs={'aniso_eq'}`.
- Lines 51: computes `offsets` using `offsets=linspace(-60e6,60e6,120)`.
- Lines 52: computes `reference_point` using `reference_point=-3.3e6`.
- Lines 55: computes `nutfrqs` using `nutfrqs=linspace(20e6,40e6,30)`.

## Implementation structure

- 2D parameter scan of a BEAM DNP experiment. <I_z> after a set
- contact time is calculated as a function of electron pulse
- amplitude and offset. Further information in:
- Calculation time: minutes (a large powder grid is needed).
- X-band magnet
- Electron and two protons
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Cartesian coordinates
- Spin temperature
- Basis set
- Spinach housekeeping
- Detection state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `kfigure()`, `offsets()`, `nutfrqs()`, `powder()`, `dnp_surf()`, `contact_curve()`, `contourf()`, `kylabel()`, `kxlabel()`, `kcolourbar()`.
