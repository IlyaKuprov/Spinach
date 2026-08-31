# examples/dnp_sol/solid_effect_field_scan_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/solid_effect_field_scan_1.m`
- Signature: `solid_effect_field_scan_1()`
- Total lines: 72

## Purpose

Magnetic field sweep DNP experiment involving a gadolinium ion, steady-state polarisation of a 15N nucleus is computed as a function of the magnetic field offset. Calculation time: minutes

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnetic field; implemented by `sys.magnet=9.4509`.
- Lines 15-16: Spin system; implemented by `sys.isotopes={'E8','15N'}`.
- Lines 18-19: Electron g-tensor; implemented by `inter.zeeman.eigs=cell(2,1)`.
- Lines 24-25: Electron ZFS tensor; implemented by `inter.coupling.eigs=cell(2,2)`.
- Lines 30-31: Coordinates (Angstrom); implemented by `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 34-35: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 39-40: Relaxation theory; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 47-48: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 51-52: Sequence parameters; implemented by `parameters.spins={'E8'}`.
- Lines 63-64: Steady state simulation; implemented by `spec=powder(spin_system,@dnp_field_scan,parameters,'esr')`.
- Lines 66-67: Plotting; implemented by `kfigure(); plot(parameters.fields,real(spec)); kgrid`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=9.4509`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'E8','15N'}`.
- Lines 19: computes `inter.zeeman.eigs` using `inter.zeeman.eigs=cell(2,1)`.
- Lines 20: computes `inter.zeeman.euler` using `inter.zeeman.euler=cell(2,1)`.
- Lines 21: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[1.9918 1.9918 1.9918]`.
- Lines 22: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[0 0 0]`.
- Lines 25: computes `inter.coupling.eigs` using `inter.coupling.eigs=cell(2,2)`.
- Lines 26: computes `inter.coupling.euler` using `inter.coupling.euler=cell(2,2)`.
- Lines 27: computes `inter.coupling.eigs{1,1}` using `inter.coupling.eigs{1,1}=570e6*[-1/3 -1/3 2*1/3]`.
- Lines 28: computes `inter.coupling.euler{1,1}` using `inter.coupling.euler{1,1}=[0 0 0]`.
- Lines 31: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 35: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 36: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 37: computes `bas.projections` using `bas.projections=[-2 -1 0 1 2]`.
- Lines 40: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 41: computes `inter.r1_rates` using `inter.r1_rates={1e4 1e1}`.
- Lines 42: computes `inter.r2_rates` using `inter.r2_rates={1e7 1e3}`.
- Lines 43: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.

## Implementation structure

- Magnetic field sweep DNP experiment involving a gadolinium ion,
- steady-state polarisation of a 15N nucleus is computed as a
- function of the magnetic field offset.
- Calculation time: minutes
- Magnetic field
- Spin system
- Electron g-tensor
- Electron ZFS tensor
- Coordinates (Angstrom)
- Basis set
- Relaxation theory
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `powder()`, `kfigure()`, `kxlabel()`, `kylabel()`.
