# examples/nmr_spen/conv_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/conv_test.m`
- Signature: `conv_test()`
- Total lines: 190

## Purpose

Convergence and accuracy test of the spatial dynamics during a pulse field gradientspin (PFG) echo sequence. The accuracy spatial diffusi- on is tested with respect to the grid size and the finite difference stancil size. The PFG spin echo pulse sequence is simulated, and then the diffusion coefficient extracted back by fitting the Stejskal-Tan- ner equation. Run time: minutes on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: Magnetic field; implemented by `sys.magnet=11.7426`.
- Lines 18-19: Isotopes; implemented by `sys.isotopes={'1H'}`.
- Lines 21-22: Chemical shift; implemented by `inter.zeeman.scalar={4.6}`.
- Lines 24-25: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 28-29: Algorithmic options; implemented by `sys.disable={'pt','krylov'}`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Aquisition parameters; implemented by `parameters.sweep=5000`.
- Lines 44-45: No relaxation; implemented by `parameters.rlx_ph={}`.
- Lines 48-49: Gradient duration; implemented by `parameters.delta_sml=0.002`.
- Lines 51-52: Diffusion delay; implemented by `parameters.delta_big=0.050`.
- Lines 54-55: Sample length; implemented by `parameters.dims=0.015`.
- Lines 57-58: Diffusion coefficient; implemented by `parameters.diff=18e-10`.
- Lines 60-61: Point stencil size list; implemented by `stencil_sizes=[3 5 7]`.
- Lines 63-64: Gradient amplitude list; implemented by `grad_amps=linspace(0,0.5,20)`.
- Lines 66-67: Grid size list; implemented by `grid_sizes=ceil(linspace(1000,10000,50))`.
- Lines 69-70: Preallocate result arrays; implemented by `D=zeros(numel(grid_sizes),numel(stencil_sizes))`.
- Lines 74-75: Grid size loop; implemented by `for k=1:numel(grid_sizes)`.
- Lines 77-78: Set grid size; implemented by `parameters.npts=grid_sizes(k)`.

### Control flow inferred from the code

- Line 75: `for` loop over `k=1:numel(grid_sizes)`.
- Line 96: `for` loop over `s=1:numel(stencil_sizes)`.
- Line 105: `parfor` loop over `n=1:numel(grad_amps)`.
- Line 156: `for` loop over `s=1:numel(stencil_sizes)`.
- Line 158: `for` loop over `k=grids_to_plot`.
- Line 165: `for` loop over `s=1:numel(stencil_sizes)`.
- Line 167: `for` loop over `k=grids_to_plot`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=11.7426`.
- Lines 19: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 22: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={4.6}`.
- Lines 25: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 26: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 29: computes `sys.disable` using `sys.disable={'pt','krylov'}`.
- Lines 30: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `parameters.sweep` using `parameters.sweep=5000`.
- Lines 38: computes `parameters.npoints` using `parameters.npoints=1024`.
- Lines 39: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 40: computes `parameters.zerofill` using `parameters.zerofill=32768`.
- Lines 41: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 42: computes `parameters.offset` using `parameters.offset=2500`.
- Lines 45: computes `parameters.rlx_ph` using `parameters.rlx_ph={}`.
- Lines 46: computes `parameters.rlx_op` using `parameters.rlx_op={}`.
- Lines 49: computes `parameters.delta_sml` using `parameters.delta_sml=0.002`.
- Lines 52: computes `parameters.delta_big` using `parameters.delta_big=0.050`.

## Implementation structure

- Convergence and accuracy test of the spatial dynamics during a pulse
- field gradientspin (PFG) echo sequence. The accuracy spatial diffusi-
- on is tested with respect to the grid size and the finite difference
- stancil size. The PFG spin echo pulse sequence is simulated, and then
- the diffusion coefficient extracted back by fitting the Stejskal-Tan-
- ner equation.
- Run time: minutes on NVidia Tesla A100, much longer on CPU
- Magnetic field
- Isotopes
- Chemical shift
- Basis set
- Algorithmic options

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `grid_sizes()`, `state()`, `report()`, `stencil_sizes()`, `tic()`, `grad_amps()`, `intensities()`, `imaging()`, `spin()`, `expfactors()`, `toc()`, `num2str()`, `kfigure()`, `scale_figure()`.
