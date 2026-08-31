# examples/esr_sol_pulsed/oop_eseem_photosystem_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/oop_eseem_photosystem_1.m`
- Signature: `oop_eseem_photosystem_1()`
- Total lines: 75

## Purpose

Powder-averaged two-pulse out-of-phase ESEEM on the [P700+,A1-] spin-correlated electron pair in Photosystem I. Time-domain si- mulation in Liouville space with averaging over a powder grid. Set to reproduce Figure 3a in The eigenvalues of P700+ g-tensor (without angles) come from The eigenvales of A-g-tensor (without angles) come from Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Magnet field; implemented by `sys.magnet=0.3249`.
- Lines 25-26: System specification; implemented by `sys.isotopes={'E','E'}`.
- Lines 36-37: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 42-43: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 46-47: Disable trajectory-level SSR algorithms; implemented by `sys.disable={'trajlevel'}`.
- Lines 49-50: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 53-54: Set the sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 63-64: Simulation; implemented by `fid=powder(spin_system,@oopeseem,parameters,'esr')`.
- Lines 66-67: Build the axis; implemented by `time_axis=(0:(parameters.npoints-1))*parameters.timestep*1e6/2`.
- Lines 69-70: Plot the results; implemented by `kfigure(); plot(time_axis,-imag(fid))`.

### Key state/data transformations

- Lines 23: computes `sys.magnet` using `sys.magnet=0.3249`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'E','E'}`.
- Lines 27: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[2.00304 2.00262 2.00232]`.
- Lines 28: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[0 0 0]`.
- Lines 29: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[2.00670 2.00560 2.00240]`.
- Lines 30: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=[0 0 0]`.
- Lines 31: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2,2)`.
- Lines 32: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=mt2hz(0.6e-3)`.
- Lines 33: computes `inter.coordinates` using `inter.coordinates={[ 0.00 0.00 0.00]`.
- Lines 37: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 38: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 39: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 40: computes `inter.damp_rate` using `inter.damp_rate=1.7e6`.
- Lines 43: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 44: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 47: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 50: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 54: computes `parameters.spins` using `parameters.spins={'E'}`.

## Implementation structure

- Powder-averaged two-pulse out-of-phase ESEEM on the [P700+,A1-]
- spin-correlated electron pair in Photosystem I. Time-domain si-
- mulation in Liouville space with averaging over a powder grid.
- Set to reproduce Figure 3a in
- The eigenvalues of P700+ g-tensor (without angles) come from
- The eigenvales of A-g-tensor (without angles) come from
- Calculation time: seconds
- Magnet field
- System specification
- Relaxation theory
- Basis set
- Disable trajectory-level SSR algorithms

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `mt2hz()`, `create()`, `basis()`, `state()`, `operator()`, `powder()`, `kfigure()`, `kxlabel()`, `kylabel()`.
