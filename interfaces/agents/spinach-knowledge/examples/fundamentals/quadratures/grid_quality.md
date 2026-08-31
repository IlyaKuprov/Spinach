# examples/fundamentals/quadratures/grid_quality.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/quadratures/grid_quality.m`
- Signature: `grid_quality()`
- Total lines: 90

## Purpose

Performance analysis for the spherical and SO(3) integration grids supplied with Spinach kernel.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-11: % Two-angle REPULSION grids; implemented by `kfigure(); hold on; kgrid; box on`.
- Lines 10-11: Create a figure; implemented by `kfigure(); hold on; kgrid; box on`.
- Lines 14-15: Loop over two-angle REPULSION grids; implemented by `for n=[100 200 400 800 1600 3200 6400 12800]`.
- Lines 17-19: Load the grid; implemented by `load(['../../../kernel/grids/rep_2ang_' num2str(n) 'pts_sph.mat'], 'alphas','betas','gammas','weights')`.
- Lines 21-22: Evaluate the grid; implemented by `grid_profile=grid_test(alphas,betas,gammas,weights,0:80,'Y_lm')`.
- Lines 24-25: Plot the evaluation; implemented by `plot(grid_profile(2:end)); axis tight`.
- Lines 28-29: Residual cosmetics; implemented by `kxlabel('spherical rank'); kylabel('integration error')`.
- Lines 34-37: % Three-angle REPULSION grids; implemented by `kfigure(); hold on; kgrid; box on`.
- Lines 40-41: Loop over three-angle REPULSION grids; implemented by `for n=[100 200 400 800 1600 3200 6400 12800]`.
- Lines 43-45: Load the grid; implemented by `load(['../../../kernel/grids/rep_3ang_' num2str(n) 'pts.mat'], 'alphas','betas','gammas','weights')`.
- Lines 47-48: Evaluate the grid; implemented by `grid_profile=grid_test(alphas,betas,gammas,weights,0:30,'D_lmn')`.
- Lines 60-63: % Two-angle Lebedev grids; implemented by `kfigure(); hold on; kgrid; box on`.
- Lines 66-67: Loop over two-angle Lebedev grids; implemented by `for n=[5 17 29 41 53]`.
- Lines 69-71: Load the grid; implemented by `load(['../../../kernel/grids/leb_2ang_rank_' num2str(n) '.mat'], 'alphas','betas','gammas','weights')`.
- Lines 73-74: Set even ranks; implemented by `ranks=2:2:100`.
- Lines 76-77: Evaluate the grid; implemented by `grid_profile=grid_test(alphas,betas,gammas,weights,ranks,'Y_lm')`.
- Lines 79-80: Plot the evaluation; implemented by `plot(ranks(2:end),grid_profile(2:end)); axis tight; drawnow`.

### Control flow inferred from the code

- Line 15: `for` loop over `n=[100 200 400 800 1600 3200 6400 12800]`.
- Line 41: `for` loop over `n=[100 200 400 800 1600 3200 6400 12800]`.
- Line 67: `for` loop over `n=[5 17 29 41 53]`.

### Key state/data transformations

- Lines 12: computes `set(gca,'YScale','log'); legend_txt` using `set(gca,'YScale','log'); legend_txt={}`.
- Lines 22: computes `grid_profile` using `grid_profile=grid_test(alphas,betas,gammas,weights,0:80,'Y_lm')`.
- Lines 26: computes `legend_txt{end+1}` using `legend_txt{end+1}=['REPULSION ' num2str(n) ' pts']`.
- Lines 74: computes `ranks` using `ranks=2:2:100`.

## Implementation structure

- Performance analysis for the spherical and SO(3) integration
- grids supplied with Spinach kernel.
- % Two-angle REPULSION grids
- Create a figure
- Loop over two-angle REPULSION grids
- Load the grid
- Evaluate the grid
- Plot the evaluation
- Residual cosmetics
- % Three-angle REPULSION grids
- Loop over three-angle REPULSION grids
- % Two-angle Lebedev grids

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `kfigure()`, `set()`, `load()`, `num2str()`, `grid_test()`, `grid_profile()`, `kxlabel()`, `kylabel()`, `klegend()`, `ranks()`.
