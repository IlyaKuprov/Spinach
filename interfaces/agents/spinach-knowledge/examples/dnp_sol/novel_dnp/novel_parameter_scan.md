# examples/dnp_sol/novel_dnp/novel_parameter_scan.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/novel_dnp/novel_parameter_scan.m`
- Signature: `novel_parameter_scan()`
- Total lines: 104

## Purpose

2D parameter scan of a NOVEL DNP experiment. <I_z> on 1H after a 0.25 us contact time is calculated as a function of electron pul- se amplitude and offset. Further information in: Calculation time: minutes

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
- Lines 54-55: Offset range, Hz; implemented by `offsets=linspace(-35e6,35e6,71)`.
- Lines 58-59: Nutation frequency range, Hz; implemented by `nutfrqs=linspace(1e6,30e6,30)`.
- Lines 61-62: Hush the output; implemented by `spin_system.sys.output='hush'`.
- Lines 64-65: Generate axis ticks; implemented by `[X,Y]=meshgrid(offsets/1e6,nutfrqs/1e6)`.
- Lines 67-68: Get the figure going; implemented by `kfigure(); dnp_surf=zeros(numel(nutfrqs),numel(offsets))`.
- Lines 70-71: Loop over powers; implemented by `for n=1:numel(nutfrqs)`.
- Lines 73-74: Set nutation frequency; implemented by `parameters.irr_powers=nutfrqs(n)`.
- Lines 76-77: Get 90-degree pulse duration; implemented by `parameters.pulse_dur=1/(4*parameters.irr_powers)`.

### Control flow inferred from the code

- Line 71: `for` loop over `n=1:numel(nutfrqs)`.
- Line 80: `parfor` loop over `m=1:numel(offsets)`.

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
- Lines 51: computes `parameters.needs` using `parameters.needs={'aniso_eq'}`.
- Lines 52: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.
- Lines 55: computes `offsets` using `offsets=linspace(-35e6,35e6,71)`.

## Implementation structure

- 2D parameter scan of a NOVEL DNP experiment. <I_z> on 1H after a
- 0.25 us contact time is calculated as a function of electron pul-
- se amplitude and offset. Further information in:
- Calculation time: minutes
- X-band magnet
- Electron and two protons
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Cartesian coordinates
- Spin temperature
- Hush the output
- Basis set
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `kfigure()`, `nutfrqs()`, `offsets()`, `powder()`, `dnp_surf()`, `contact_curve()`, `contourf()`, `kylabel()`, `kxlabel()`, `kcolourbar()`.
