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
