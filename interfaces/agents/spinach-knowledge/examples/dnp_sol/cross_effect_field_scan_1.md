# examples/dnp_sol/cross_effect_field_scan_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/cross_effect_field_scan_1.m`
- Signature: `cross_effect_field_scan_1()`
- Total lines: 82

## Purpose

Magnetic field sweep cross effect DNP experiment -steady-state proton magnetisation under microwave iradiation as a function of the applied magnetic field. A powder average calculation. Calculation time: hours

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Magnetic field; implemented by `sys.magnet=18.78`.
- Lines 14-15: Spin system; implemented by `sys.isotopes={'E','E','14N','1H'}`.
- Lines 17-18: Electron g-tensors; implemented by `inter.zeeman.eigs=cell(4,1)`.
- Lines 25-26: 14N quadrupolar tensor; implemented by `inter.coupling.eigs=cell(4,4)`.
- Lines 31-32: Coordinates (Angstrom); implemented by `inter.coordinates={[ 0.00 0.00 0.00]`.
- Lines 37-38: Hyperfine couplings; implemented by `inter.coupling.eigs{1,3}=[17.4e6 17.6e6 102e6]`.
- Lines 41-42: Exchange coupling; implemented by `inter.coupling.scalar=cell(4,4)`.
- Lines 45-46: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 49-50: Relaxation theory; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 57-58: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 61-62: Sequence parameters; implemented by `parameters.mw_pwr=10e6`.
- Lines 73-74: Steady state simulation; implemented by `answer=powder(spin_system,@dnp_field_scan,parameters,'esr')`.
- Lines 76-77: Plotting; implemented by `kfigure(); plot(parameters.fields,real(answer)); kgrid`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=18.78`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'E','E','14N','1H'}`.
- Lines 18: computes `inter.zeeman.eigs` using `inter.zeeman.eigs=cell(4,1)`.
- Lines 19: computes `inter.zeeman.euler` using `inter.zeeman.euler=cell(4,1)`.
- Lines 20: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[2.0085 2.00605 2.00215]`.
- Lines 21: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[pi/2 pi/3 pi/4]`.
- Lines 22: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[2.00319 2.00319 2.00258]`.
- Lines 23: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=[pi/5 pi/6 pi/7]`.
- Lines 26: computes `inter.coupling.eigs` using `inter.coupling.eigs=cell(4,4)`.
- Lines 27: computes `inter.coupling.euler` using `inter.coupling.euler=cell(4,4)`.
- Lines 28: computes `inter.coupling.eigs{3,3}` using `inter.coupling.eigs{3,3}=[-1e6 -1e6 2e6]`.
- Lines 29: computes `inter.coupling.euler{3,3}` using `inter.coupling.euler{3,3}=[0 0 0]`.
- Lines 32: computes `inter.coordinates` using `inter.coordinates={[ 0.00 0.00 0.00]`.
- Lines 38: computes `inter.coupling.eigs{1,3}` using `inter.coupling.eigs{1,3}=[17.4e6 17.6e6 102e6]`.
- Lines 39: computes `inter.coupling.euler{1,3}` using `inter.coupling.euler{1,3}=[0 0 0]`.
- Lines 42: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(4,4)`.
- Lines 43: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=2*(-73e6)`.
- Lines 46: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.

## Implementation structure

- Magnetic field sweep cross effect DNP experiment -steady-state proton
- magnetisation under microwave iradiation as a function of the applied
- magnetic field. A powder average calculation.
- Calculation time: hours
- Magnetic field
- Spin system
- Electron g-tensors
- 14N quadrupolar tensor
- Coordinates (Angstrom)
- Hyperfine couplings
- Exchange coupling
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `powder()`, `kfigure()`, `kxlabel()`, `kylabel()`.
