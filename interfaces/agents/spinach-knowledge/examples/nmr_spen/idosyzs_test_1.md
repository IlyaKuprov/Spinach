# examples/nmr_spen/idosyzs_test_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/idosyzs_test_1.m`
- Signature: `idosyzs_test_1()`
- Total lines: 104

## Purpose

Diffusion attenuation during soft pulses in a simplified model sequence of the Zangger-Sterk pure shift iDOSY with fitting using a modified ver- sion of the Stejskal Tanner equation, as described in: Calculation time: seconds on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: Magnetic field; implemented by `sys.magnet=11.7426`.
- Lines 18-19: Isotopes; implemented by `sys.isotopes={'1H'}`.
- Lines 21-22: Chemical shift; implemented by `inter.zeeman.scalar={4.6}`.
- Lines 24-25: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 28-29: Algorithmic options; implemented by `sys.disable={'pt'}`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Sample geometry; implemented by `parameters.dims=0.015`.
- Lines 41-42: Diffusion coefficient; implemented by `parameters.diff=18e-10`.
- Lines 44-45: Relaxation phantom; implemented by `parameters.rlx_ph={}`.
- Lines 48-49: White margins on the initial condition; implemented by `parameters.rho0_ph={[zeros(1000,1); ones(2000,1); zeros(1000,1)]}`.
- Lines 54-55: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 65-66: Gradient amplitude list; implemented by `grad_amps=linspace(0.01,0.40,20)`.
- Lines 68-69: Compute a series of spectra with different gradient amplitudes; implemented by `parfor n=1:numel(grad_amps)`.
- Lines 71-72: Set gradient amplitude; implemented by `localpar=parameters; localpar.g_amp=grad_amps(n)`.
- Lines 74-75: Run the simulation; implemented by `inten(n)=imaging(spin_system,@idosyzs,localpar)`.
- Lines 79-80: Normalise and gather intensities; implemented by `inten=gather(inten/inten(1))`.
- Lines 82-84: Stejskal-Tanner fitting function; implemented by `expfactor=(spin(parameters.spins{1})*parameters.delta_sml).^2* (parameters.delta_big-(parameters.delta_sml)/3)*1e-10`.
- Lines 87-88: Run Stejskal-Tanner fitting; implemented by `opts=optimoptions(@lsqcurvefit,'Display','iter')`.

### Control flow inferred from the code

- Line 69: `parfor` loop over `n=1:numel(grad_amps)`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=11.7426`.
- Lines 19: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 22: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={4.6}`.
- Lines 25: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 26: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 29: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 30: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `parameters.dims` using `parameters.dims=0.015`.
- Lines 38: computes `parameters.npts` using `parameters.npts=4000`.
- Lines 39: computes `parameters.deriv` using `parameters.deriv={'period',7}`.
- Lines 42: computes `parameters.diff` using `parameters.diff=18e-10`.
- Lines 45: computes `parameters.rlx_ph` using `parameters.rlx_ph={}`.
- Lines 46: computes `parameters.rlx_op` using `parameters.rlx_op={}`.
- Lines 49: computes `parameters.rho0_ph` using `parameters.rho0_ph={[zeros(1000,1); ones(2000,1); zeros(1000,1)]}`.
- Lines 50: computes `parameters.rho0_st` using `parameters.rho0_st={state(spin_system,'Lz','1H')}`.
- Lines 51: computes `parameters.coil_ph` using `parameters.coil_ph={ones(parameters.npts,1)}`.
- Lines 52: computes `parameters.coil_st` using `parameters.coil_st={state(spin_system,'L+','1H')}`.

## Implementation structure

- Diffusion attenuation during soft pulses in a simplified model sequence
- of the Zangger-Sterk pure shift iDOSY with fitting using a modified ver-
- sion of the Stejskal Tanner equation, as described in:
- Calculation time: seconds on NVidia Tesla A100, much longer on CPU
- Magnetic field
- Isotopes
- Chemical shift
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sample geometry
- Diffusion coefficient

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `grad_amps()`, `inten()`, `imaging()`, `gather()`, `spin()`, `optimoptions()`, `lsqcurvefit()`, `kfigure()`, `num2str()`.
